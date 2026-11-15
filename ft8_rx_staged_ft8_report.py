import numpy as np
import pmt
from gnuradio import gr

FS = 12000                      # decoder sample rate
NSPS = 1920                     # samples per FT8 symbol (0.16 s)
NN = 79                         # symbols per transmission
COSTAS = np.array([3, 1, 4, 0, 6, 5, 2])
GRAY = np.array([0, 1, 3, 2, 5, 6, 4, 7])          # 3 data bits -> tone
DATA_SYMS = np.r_[7:36, 43:72]
SYNC_SYMS = np.r_[0:7, 36:43, 72:79]


# ---- message helpers: every stage passes one PMT dict per slot and pass ----
def to_msg(d):
    """Python dict -> PMT dict. numpy arrays become uniform vectors (fast path)."""
    m = pmt.make_dict()
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            a = v.ravel()
            if a.dtype == np.uint8:
                p = pmt.init_u8vector(len(a), a.tolist())
            elif a.dtype.kind == 'i':
                p = pmt.init_s32vector(len(a), a.astype(np.int32).tolist())
            elif a.dtype.kind == 'c':
                p = pmt.init_c32vector(len(a), a.astype(np.complex64).tolist())
            elif a.dtype == np.float64:
                p = pmt.init_f64vector(len(a), a.tolist())
            else:
                p = pmt.init_f32vector(len(a), a.astype(np.float32).tolist())
        else:
            p = pmt.to_pmt(v)
        m = pmt.dict_add(m, pmt.intern(k), p)
    return m


def from_msg(m):
    return pmt.to_python(m)


# ---------------------------------------------------------------- 77-bit message unpack
_A37 = ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_A36 = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_A27 = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_A38 = ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ/'
_A42 = ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+-./?'
_NTOKENS = 2063592
_MAX22 = 4194304
_MAXGRID4 = 32400
_CSEC = ('AB AK AL AR AZ BC CO CT DE EB EMA ENY EPA EWA GA GTA IA ID IL IN KS KY LA LAX MAR MB MDC '
         'ME MI MN MO MS MT NC ND NE NFL NH NL NLI NM NNJ NNY NT NTX NV OH OK ONE ONN ONS OR ORG PAC '
         'PR QC RI SB SC SCV SD SDG SF SFL SJV SK SNJ STX SV TN UT VA VI VT WCF WI WMA WNY WPA WTX '
         'WV WWA WY DX').split()
_CMULT = ('AL AK AZ AR CA CO CT DE FL GA HI ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ '
          'NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY NB NS QC ON MB SK AB BC NWT NF '
          'LB NU YT PEI DC').split()


class CallHash:
    """Callsigns heard so far, keyed by their 22/12/10-bit hashes (for <...> fields)."""

    def __init__(self, max_entries=2000):
        self.h22, self.h12, self.h10 = {}, {}, {}
        self.max_entries = max_entries

    @staticmethod
    def hash22(call):
        n = 0
        c = (call + ' ' * 11)[:11]
        for ch in c:
            j = _A38.find(ch)
            if j < 0:
                return None
            n = 38 * n + j
        return (((47055833459 * n) & 0xFFFFFFFFFFFFFFFF) >> 42) & 0x3FFFFF

    def save(self, call):
        call = call.strip()
        if len(call) < 3 or call.startswith('<'):
            return
        h = self.hash22(call)
        if h is None:
            return
        if len(self.h22) > self.max_entries:
            self.h22.clear(); self.h12.clear(); self.h10.clear()
        self.h22[h] = call
        self.h12[h >> 10] = call
        self.h10[h >> 12] = call

    def lookup(self, h, nbits):
        d = {22: self.h22, 12: self.h12, 10: self.h10}[nbits]
        c = d.get(h)
        return '<%s>' % c if c else '<...>'


def _bits_to_int(b):
    v = 0
    for x in b:
        v = (v << 1) | int(x)
    return v


def _unpack28(n28, ip, i3, hashes):
    if n28 < _NTOKENS:
        if n28 == 0:
            return 'DE'
        if n28 == 1:
            return 'QRZ'
        if n28 == 2:
            return 'CQ'
        if n28 <= 1002:
            return 'CQ %03d' % (n28 - 3)
        if n28 <= 532443:
            n = n28 - 1003
            s = ''
            for _ in range(4):
                s = _A27[n % 27] + s
                n //= 27
            return 'CQ ' + s.strip()
        return None
    n28 -= _NTOKENS
    if n28 < _MAX22:
        return hashes.lookup(n28, 22)
    n = n28 - _MAX22
    c = [''] * 6
    c[5] = _A27[n % 27]; n //= 27
    c[4] = _A27[n % 27]; n //= 27
    c[3] = _A27[n % 27]; n //= 27
    c[2] = '0123456789'[n % 10]; n //= 10
    c[1] = _A36[n % 36]; n //= 36
    if n >= 37:
        return None
    c[0] = _A37[n]
    call = ''.join(c)
    if call.startswith('3D0') and call[3] != ' ':
        call = '3DA0' + call[3:].strip()
    elif call[0] == 'Q' and call[1].isalpha():
        call = '3X' + call[1:].strip()
    call = call.strip()
    if len(call) < 3 or ' ' in call:
        return None
    hashes.save(call)
    if ip:
        if i3 == 1:
            call += '/R'
        elif i3 == 2:
            call += '/P'
        else:
            return None
        hashes.save(call)
    return call


def _grid_or_report(igrid4, ir):
    if igrid4 <= _MAXGRID4:
        n = igrid4
        d3 = n % 10; n //= 10
        d2 = n % 10; n //= 10
        l1 = n % 18; n //= 18
        l0 = n % 18
        g = chr(65 + l0) + chr(65 + l1) + str(d2) + str(d3)
        return ('R ' + g) if ir else g
    irpt = igrid4 - _MAXGRID4
    if irpt == 1:
        return ''
    if irpt == 2:
        return 'RRR'
    if irpt == 3:
        return 'RR73'
    if irpt == 4:
        return '73'
    r = '%+03d' % (irpt - 35)
    return ('R' + r) if ir else r


def unpack77(b, hashes):
    """b: 77 bits (array of 0/1). Returns message text or None (unsupported/invalid)."""
    b = [int(x) for x in b]
    i3 = _bits_to_int(b[74:77])
    n3 = _bits_to_int(b[71:74])
    if i3 in (1, 2):
        c1 = _unpack28(_bits_to_int(b[0:28]), b[28], i3, hashes)
        c2 = _unpack28(_bits_to_int(b[29:57]), b[57], i3, hashes)
        if c1 is None or c2 is None:
            return None
        ex = _grid_or_report(_bits_to_int(b[59:74]), b[58])
        if c1.startswith('CQ') and ex.startswith('R'):
            return None                      # implausible: CQ with R/report
        return ' '.join(x for x in (c1, c2, ex) if x)
    if i3 == 4:
        n12 = _bits_to_int(b[0:12])
        n58 = _bits_to_int(b[12:70])
        iflip, nrpt, icq = b[70], _bits_to_int(b[71:73]), b[73]
        s = ''
        for _ in range(11):
            s = _A38[n58 % 38] + s
            n58 //= 38
        full = s.strip()
        if len(full) < 3 or ' ' in full:
            return None
        hashes.save(full)
        hc = hashes.lookup(n12, 12)
        if icq:
            return 'CQ ' + full
        c1, c2 = (full, hc) if iflip else (hc, full)
        ex = ('', 'RRR', 'RR73', '73')[nrpt]
        return ' '.join(x for x in (c1, c2, ex) if x)
    if i3 == 3:
        tu = b[0]
        c1 = _unpack28(_bits_to_int(b[1:29]), 0, 3, hashes)
        c2 = _unpack28(_bits_to_int(b[29:57]), 0, 3, hashes)
        if c1 is None or c2 is None:
            return None
        r = b[57]
        rst = '5%d9' % (_bits_to_int(b[58:61]) + 2)
        s13 = _bits_to_int(b[61:74])
        if s13 > 8000:
            if s13 - 8001 >= len(_CMULT):
                return None
            ex = _CMULT[s13 - 8001]
        else:
            ex = '%04d' % s13
        parts = (['TU;'] if tu else []) + [c1, c2] + (['R'] if r else []) + [rst, ex]
        return ' '.join(parts)
    if i3 == 5:
        h12 = _bits_to_int(b[0:12])
        h22 = _bits_to_int(b[12:34])
        r = b[34]
        rst = '5%d9' % (_bits_to_int(b[35:38]) + 2)
        serial = _bits_to_int(b[38:49])
        n = _bits_to_int(b[49:74])
        d = [18 * 10 * 10 * 24 * 24, 10 * 10 * 24 * 24, 10 * 24 * 24, 24 * 24, 24, 1]
        j = []
        for dv in d:
            j.append(n // dv)
            n -= j[-1] * dv
        if j[0] > 17 or j[1] > 17 or j[4] > 23 or j[5] > 23:
            return None
        g6 = chr(65 + j[0]) + chr(65 + j[1]) + str(j[2]) + str(j[3]) + chr(65 + j[4]) + chr(65 + j[5])
        return '%s %s %s%s%04d %s' % (hashes.lookup(h12, 12), hashes.lookup(h22, 22),
                                       'R ' if r else '', rst, serial, g6)
    if i3 == 0:
        if n3 == 0:                          # free text
            n = _bits_to_int(b[0:71])
            s = ''
            for _ in range(13):
                s = _A42[n % 42] + s
                n //= 42
            s = s.strip()
            return s if s else None
        if n3 == 1:                          # DXpedition
            c1 = _unpack28(_bits_to_int(b[0:28]), 0, 0, hashes)
            c2 = _unpack28(_bits_to_int(b[28:56]), 0, 0, hashes)
            if c1 is None or c2 is None:
                return None
            h10 = _bits_to_int(b[56:66])
            rpt = 2 * _bits_to_int(b[66:71]) - 30
            return '%s RR73; %s %s %+03d' % (c1, c2, hashes.lookup(h10, 10), rpt)
        if n3 in (3, 4):                     # ARRL Field Day
            c1 = _unpack28(_bits_to_int(b[0:28]), 0, 0, hashes)
            c2 = _unpack28(_bits_to_int(b[28:56]), 0, 0, hashes)
            if c1 is None or c2 is None:
                return None
            r = b[56]
            ntx = _bits_to_int(b[57:61]) + (1 if n3 == 3 else 17)
            cls = chr(65 + _bits_to_int(b[61:64]))
            isec = _bits_to_int(b[64:71])
            if isec < 1 or isec > len(_CSEC):
                return None
            return ' '.join([c1, c2] + (['R'] if r else []) + ['%d%s' % (ntx, cls), _CSEC[isec - 1]])
        if n3 == 5:                          # telemetry
            v = _bits_to_int(b[0:71])
            return ('%018X' % v).lstrip('0') or '0'
    return None


# ================================================================ Stage 5 block
import sys
import time


def tones_from_cw(cw):
    """174-bit codeword -> 79 channel tones (Costas + Gray-mapped data)."""
    cw = np.asarray(cw, dtype=int)
    t = np.zeros(NN, dtype=int)
    t[SYNC_SYMS] = np.tile(COSTAS, 3)
    t[DATA_SYMS] = GRAY[cw[0::3] * 4 + cw[1::3] * 2 + cw[2::3]]
    return t


class blk(gr.basic_block):
    """FT8 Unpack + Report: codewords -> messages, table, and the subtraction loop.

    Unpacks each new 77-bit payload (duplicates from earlier passes are ignored),
    measures SNR (2500 Hz reference), DT and frequency, and keeps the slot's list.
    If more passes are allowed and this pass found something new, the new decodes go
    to FT8 Subtract ('new'), which sends the residual back to the Sync Search.
    Otherwise the slot is finished: WSJT-X style table on the console and in the
    QT window (Python Snippet calls make_widget()), one PMT dict per decode on 'decodes'.
    """

    def __init__(self, passes=3, console=True):
        gr.basic_block.__init__(self, name='FT8 Unpack + Report', in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern('cw'))
        self.message_port_register_out(pmt.intern('new'))
        self.message_port_register_out(pmt.intern('decodes'))
        self.set_msg_handler(pmt.intern('cw'), self.handle)
        self.passes = max(1, int(passes))
        self.console = bool(console)
        self.hashes = CallHash()
        self.slots = {}
        self._bridge = None
        self._header_done = False

    def make_widget(self):
        """Create the decode table. Call from the GUI thread (Python Snippet, Main - After Init)."""
        from PyQt5 import QtCore, QtGui, QtWidgets

        class _Bridge(QtCore.QObject):
            slot_done = QtCore.pyqtSignal(object)

        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        top = QtWidgets.QHBoxLayout()
        status = QtWidgets.QLabel('Waiting for the first FT8 slot ...')
        clear = QtWidgets.QPushButton('Clear')
        top.addWidget(status, 1)
        top.addWidget(clear)
        lay.addLayout(top)
        tbl = QtWidgets.QTableWidget(0, 5)
        tbl.setHorizontalHeaderLabels(['UTC', 'dB', 'DT', 'Freq', 'Message'])
        tbl.verticalHeader().setVisible(False)
        tbl.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        tbl.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        font = QtGui.QFont('Consolas')
        font.setStyleHint(QtGui.QFont.Monospace)
        tbl.setFont(font)
        hh = tbl.horizontalHeader()
        fm = QtGui.QFontMetrics(font)
        for c, sample in enumerate(('000000', '-24', '-1.5', '0000')):
            hh.setSectionResizeMode(c, QtWidgets.QHeaderView.Interactive)
            tbl.setColumnWidth(c, fm.width(sample) + 18)
        hh.setSectionResizeMode(4, QtWidgets.QHeaderView.Stretch)
        tbl.verticalHeader().setDefaultSectionSize(fm.height() + 4)
        lay.addWidget(tbl)
        w.setMinimumWidth(460)
        clear.clicked.connect(lambda: tbl.setRowCount(0))
        cq_bg = QtGui.QBrush(QtGui.QColor(200, 255, 200))
        sep_bg = QtGui.QBrush(QtGui.QColor(225, 225, 225))

        def on_slot(info):
            sb = tbl.verticalScrollBar()
            at_bottom = sb.value() >= sb.maximum() - 2
            r = tbl.rowCount()
            tbl.insertRow(r)
            it = QtWidgets.QTableWidgetItem('---- %s  %d decodes  (%.1f s)' % (info['label'], len(info['rows']), info['secs']))
            it.setBackground(sep_bg)
            tbl.setItem(r, 0, it)
            tbl.setSpan(r, 0, 1, 5)
            for d in info['rows']:
                r = tbl.rowCount()
                tbl.insertRow(r)
                vals = (info['utc'], '%d' % d['snr'], '%.1f' % d['dt'], '%d' % round(d['freq']), d['msg'])
                for c, v in enumerate(vals):
                    item = QtWidgets.QTableWidgetItem(v)
                    if c < 4:
                        item.setTextAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
                    if d['msg'].startswith('CQ '):
                        item.setBackground(cq_bg)
                    tbl.setItem(r, c, item)
            excess = tbl.rowCount() - 3000
            for _ in range(max(0, excess)):
                tbl.removeRow(0)
            if at_bottom:
                tbl.scrollToBottom()
            status.setText(info['status'])

        self._bridge = _Bridge()
        self._bridge.slot_done.connect(on_slot)
        return w

    def handle(self, m):
        d = from_msg(m)
        key = d['slot']
        st = self.slots.setdefault(key, {'seen': set(), 'rows': []})
        k = len(d['f1'])
        cw = np.asarray(d['cw']).reshape(k, 174) if k else np.zeros((0, 174), np.uint8)
        cs = np.asarray(d['cs']).reshape(k, NN, 8) if k else None
        base = np.asarray(d['base'])
        new = []
        for c in range(k):
            a77 = cw[c, :77].astype(np.uint8)
            kk = a77.tobytes()
            if kk in st['seen']:
                continue
            msg = unpack77(a77, self.hashes)
            if msg is None:
                continue
            st['seen'].add(kk)
            f1 = float(d['f1'][c])
            tones = tones_from_cw(cw[c])
            sig_p = (np.abs(cs[c][np.arange(NN), tones]) ** 2).mean()
            i0 = max(0, int((f1 - 300.0) / 3.125))
            i1 = min(len(base), max(i0 + 2, int((f1 + 350.0) / 3.125)))
            noise_cs = 32.0 * (np.percentile(base[i0:i1], 20) / 0.513) / NSPS
            arg = sig_p / max(noise_cs, 1e-30) - 1.0
            snr = 10.0 * np.log10(max(arg, 1e-3)) - 28.3
            new.append(dict(snr=int(round(max(snr, -24.0))), dt=d['ib'][c] / 200.0 - 0.5,
                            freq=f1, msg=msg, tones=tones))
        st['rows'] += new
        if new and d['pass'] + 1 < self.passes:
            self.message_port_pub(pmt.intern('new'), to_msg(dict(
                slot=d['slot'], label=d['label'], pretty=d['pretty'], **{'pass': d['pass']}, t0=d['t0'],
                base=base, tones=np.array([r['tones'] for r in new], dtype=np.uint8).reshape(-1),
                f1=np.array([r['freq'] for r in new], dtype=np.float64),
                tstart=np.array([r['dt'] + 0.5 for r in new], dtype=np.float64))))
            return
        self._finish(d, self.slots.pop(key)['rows'])

    def _finish(self, d, rows):
        rows.sort(key=lambda r: r['freq'])
        secs = time.time() - d['t0']
        utc, pretty = d['label'], d['pretty']
        if self.console:
            out = []
            if not self._header_done:
                out.append('UTC     dB   DT  Freq  Message')
                self._header_done = True
            out.append('------ %s : %d decodes (%.1f s, %d pass%s) ------'
                       % (pretty, len(rows), secs, d['pass'] + 1, '' if d['pass'] == 0 else 'es'))
            for r in rows:
                out.append('%s %3d %4.1f %5d ~  %s' % (utc, r['snr'], r['dt'], round(r['freq']), r['msg']))
            sys.stdout.write('\n'.join(out) + '\n')
            sys.stdout.flush()
        for r in rows:
            meta = pmt.make_dict()
            meta = pmt.dict_add(meta, pmt.intern('utc'), pmt.intern(utc))
            meta = pmt.dict_add(meta, pmt.intern('snr'), pmt.from_long(int(r['snr'])))
            meta = pmt.dict_add(meta, pmt.intern('dt'), pmt.from_double(float(r['dt'])))
            meta = pmt.dict_add(meta, pmt.intern('freq'), pmt.from_double(float(r['freq'])))
            meta = pmt.dict_add(meta, pmt.intern('msg'), pmt.intern(r['msg']))
            self.message_port_pub(pmt.intern('decodes'), meta)
        if self._bridge is not None:
            status = 'Last slot %s: %d decodes, %.1f s, %d pass%s' % (
                pretty, len(rows), secs, d['pass'] + 1, '' if d['pass'] == 0 else 'es')
            self._bridge.slot_done.emit(dict(label=pretty, utc=utc, rows=rows, secs=secs, status=status))
