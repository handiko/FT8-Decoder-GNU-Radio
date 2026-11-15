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


# ================================================================ Stage 1: slot clock
import sys
import time
import datetime

_SLOT = 15.0


class blk(gr.sync_block):
    """FT8 Slot Clock: 12 kHz audio in -> one message per 15 s slot (the slot's audio).

    pc_clock=True : slots follow the PC clock (UTC :00/:15/:30/:45); live sound card.
    pc_clock=False: slots counted from the first sample (WSJT-X .wav test files).
    decode_at     : seconds into the slot when the slot is sent (live mode).
    clock_offset  : seconds added to the PC clock (sound-card latency), usually 0.
    Output 'slot': {slot, label, pretty, pass=0, t0, audio[180000]}
    """

    def __init__(self, samp_rate=12000, decode_at=14.4, pc_clock=True, clock_offset=0.0):
        gr.sync_block.__init__(self, name='FT8 Slot Clock', in_sig=[np.float32], out_sig=None)
        self.message_port_register_out(pmt.intern('slot'))
        if int(samp_rate) != FS:
            raise ValueError('FT8 Slot Clock needs 12000 S/s input')
        self.decode_at = float(decode_at)
        self.pc_clock = bool(pc_clock)
        self.clock_offset = float(clock_offset)
        self.buf = np.zeros(40 * FS, dtype=np.float32)        # ring buffer
        self.n_total = 0
        self.n_first = None
        self.t_anchor = None                                    # PC time of sample n_anchor
        self.n_anchor = 0
        self.err_min = np.inf
        self.last_slot = None

    def work(self, input_items, output_items):
        x = input_items[0]
        n = len(x)
        now = time.time() + self.clock_offset
        L = len(self.buf)
        i = self.n_total % L
        if n >= L:
            self.buf[:] = x[-L:]
        else:
            k = min(n, L - i)
            self.buf[i:i + k] = x[:k]
            self.buf[:n - k] = x[k:]
        if self.n_first is None:
            self.n_first = self.n_total
        self.n_total += n

        if self.pc_clock:
            if self.t_anchor is None:
                self.t_anchor, self.n_anchor = now, self.n_total
            err = now - self._t_of(self.n_total)      # >= 0 apart from drift
            if abs(err) > 1.0:                         # lost samples / clock step: re-anchor
                self.t_anchor, self.n_anchor = now, self.n_total
                self.err_min = np.inf
                err = 0.0
            self.err_min = min(self.err_min, err)
            ts = np.floor((self._t_of(self.n_total) - self.decode_at) / _SLOT) * _SLOT
            if self.last_slot is None:
                self.last_slot = ts                    # never decode a slot we joined late
            elif ts > self.last_slot:
                self.last_slot = ts
                if np.isfinite(self.err_min):           # follow sound-card clock drift
                    self.t_anchor += self.err_min
                    self.err_min = np.inf
                n0 = self.n_anchor + int(round((ts - self.t_anchor) * FS))
                dtu = datetime.datetime.fromtimestamp(ts, datetime.timezone.utc)
                self._send(n0, ts, dtu.strftime('%H%M%S'), dtu.strftime('%Y-%m-%d %H:%M:%S UTC'))
        else:
            k = self.n_total // (15 * FS)
            if self.last_slot is None:
                self.last_slot = -1
            if k - 1 > self.last_slot:
                self.last_slot = k - 1
                ts = self.last_slot * _SLOT
                lab = '%02d%02d%02d' % (ts // 3600 % 24, ts // 60 % 60, ts % 60)
                self._send(self.n_first + self.last_slot * 15 * FS, ts, lab, 'slot %s (file time)' % lab)
        return n

    def _t_of(self, nn):
        return self.t_anchor + (nn - self.n_anchor) / float(FS)

    def _send(self, n0, ts, label, pretty):
        L = len(self.buf)
        seg = np.zeros(15 * FS, dtype=np.float32)
        lo = max(n0, self.n_total - L, self.n_first)
        hi = min(n0 + 15 * FS, self.n_total)
        if hi - lo < 13 * FS:                           # less than 13 s of this slot: skip
            return
        seg[lo - n0:hi - n0] = self.buf[np.arange(lo, hi) % L]
        self.message_port_pub(pmt.intern('slot'), to_msg(dict(
            slot=float(ts), label=label, pretty=pretty, **{'pass': 0}, t0=time.time(), audio=seg)))
