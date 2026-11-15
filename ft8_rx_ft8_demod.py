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


# ---------------------------------------------------------------- DSP
_FS2 = 200
_NFFT_BIG = 192000
_NFFT2 = 3200
_TAPER = 0.5 * (1.0 + np.cos(np.arange(101) * np.pi / 100.0))
_M32 = np.arange(32)
_CSYNC = np.exp(2j * np.pi * np.outer(COSTAS, _M32) / 32.0)
_SYNC_OFF = np.concatenate([32 * (np.arange(7) + 36 * g) for g in range(3)])
_SYNC_TONE = np.tile(COSTAS, 3)


class Downsampler:
    """One big FFT of the slot; per call: band around f0 -> 200 S/s complex at 0 Hz."""

    def __init__(self, x):
        self.cx = np.fft.rfft(x, _NFFT_BIG)
        self.df = FS / _NFFT_BIG

    def __call__(self, f0):
        df = self.df
        i0 = int(round(f0 / df))
        it = min(int(round((f0 + 8.5 * 6.25) / df)), _NFFT_BIG // 2)
        ib = max(1, int(round((f0 - 1.5 * 6.25) / df)))
        k = it - ib + 1
        c1 = np.zeros(_NFFT2, dtype=complex)
        c1[:k] = self.cx[ib:it + 1]
        c1[:101] *= _TAPER[::-1]
        c1[k - 101:k] *= _TAPER
        c1 = np.roll(c1, -(i0 - ib))
        return np.fft.ifft(c1) * (_NFFT2 / np.sqrt(_NFFT_BIG * _NFFT2))


def _sync_power(cdp, starts, ref):
    """Coherent Costas correlation power for each start index (into padded cd)."""
    idx = starts[:, None, None] + _SYNC_OFF[None, :, None] + _M32[None, None, :]
    seg = cdp[idx]                                          # (P,21,32)
    z = (seg * ref[None]).sum(axis=2)
    return (z.real ** 2 + z.imag ** 2).sum(axis=1)


_PADL = 400                                                 # padding (samples @200 S/s)


def fine_sync(ds, f0, xdt):
    """WSJT-X style refinement. Returns (f1, i_start, cd) with cd padded by _PADL."""
    cd = ds(f0)
    cdp = np.concatenate([np.zeros(_PADL, complex), cd, np.zeros(_PADL + 2600, complex)])
    i0 = int(round((xdt + 0.5) * _FS2))
    ref0 = np.conj(np.tile(_CSYNC, (3, 1)))                # (21,32)
    starts = np.arange(i0 - 10, i0 + 11) + _PADL
    starts = np.clip(starts, 0, len(cdp) - 1 - 32 * 79)
    p = _sync_power(cdp, starts, ref0)
    ibest = int(starts[np.argmax(p)]) - _PADL
    best_f, best_p = 0.0, -1.0
    tw_t = _M32 / _FS2
    for delf in np.arange(-2.5, 2.51, 0.5):
        ref = ref0 * np.exp(-2j * np.pi * delf * tw_t)[None, :]
        pp = _sync_power(cdp, np.array([ibest + _PADL]), ref)[0]
        if pp > best_p:
            best_p, best_f = pp, delf
    f1 = f0 + best_f
    cd = ds(f1)
    cdp = np.concatenate([np.zeros(_PADL, complex), cd, np.zeros(_PADL + 2600, complex)])
    starts = np.arange(ibest - 4, ibest + 5) + _PADL
    starts = np.clip(starts, 0, len(cdp) - 1 - 32 * 79)
    p = _sync_power(cdp, starts, ref0)
    ibest = int(starts[np.argmax(p)]) - _PADL
    return f1, ibest, cdp


def _bit_tables():
    """For nsym=1,2,3: tone index tuples and bit masks for all 8**nsym hypotheses."""
    tabs = {}
    for nsym in (1, 2, 3):
        nt = 8 ** nsym
        i = np.arange(nt)
        digs = [(i >> (3 * (nsym - 1 - d))) & 7 for d in range(nsym)]   # per symbol 3-bit value
        tone = np.stack([GRAY[dg] for dg in digs])                       # (nsym, nt)
        nb = 3 * nsym
        one = np.stack([(i >> (nb - 1 - b)) & 1 for b in range(nb)]).astype(bool)  # (nb, nt) MSB first
        tabs[nsym] = (tone, one)
    return tabs


_BITTAB = _bit_tables()


def soft_metrics(cs):
    """cs: (79,8) complex tone amplitudes -> (4,174) normalized LLRs.
    Rows: 1-symbol, 2-symbol and 3-symbol coherent metrics, and the 1-symbol
    metric normalized per bit (same four variants WSJT-X tries)."""
    out = np.zeros((4, 174))
    for nsym in (1, 2, 3):
        tone, one = _BITTAB[nsym]
        nb = 3 * nsym
        for ihalf in range(2):
            k = np.arange(0, 29, nsym)
            ks = k + 7 + 36 * ihalf
            acc = np.zeros((len(ks), tone.shape[1]), dtype=complex)
            for d in range(nsym):
                kk = ks + d
                acc += np.where((kk < 79)[:, None], cs[np.minimum(kk, 78)][:, tone[d]], 0)
            s2 = np.abs(acc)[:, None, :]                       # (G,1,nt)
            m1 = np.where(one[None], s2, -1.0).max(axis=2)     # (G,nb)
            m0 = np.where(~one[None], s2, -1.0).max(axis=2)
            pos = (k * 3 + ihalf * 87)[:, None] + np.arange(nb)[None, :]
            ok = pos < 174
            out[nsym - 1, pos[ok]] = (m1 - m0)[ok]
            if nsym == 1:
                den = np.maximum(m1, m0)
                out[3, pos[ok]] = np.where(den > 0, (m1 - m0) / np.where(den > 0, den, 1), 0)[ok]
    sd = out.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return 2.83 * out / sd


def symbol_spectra(cdp, ibest):
    idx = ibest + _PADL + 32 * np.arange(NN)[:, None] + _M32[None, :]
    seg = cdp[idx]                                          # (79,32)
    return np.fft.fft(seg, axis=1)[:, :8]                   # (79,8)


# ================================================================ Stage 3 block
class blk(gr.basic_block):
    """FT8 Demodulator: candidates -> soft bits.

    Per candidate: band-limit, shift to 0 Hz and decimate to 200 S/s (one FFT of the
    whole slot, one small IFFT per candidate), fine time (+-50 ms) and frequency
    (+-2.5 Hz) sync on the Costas arrays, 79 x 8 tone amplitudes, and four soft-bit
    metrics (1-, 2-, 3-symbol coherent and normalized). Candidates whose sync tones are
    mostly wrong (<= 6 of 21) are dropped.
    Output 'llr'  : slot info + f1[], ib[], llr[N*4*174], cs[N*79*8], base[]
    Output 'audio': slot, pass, audio (for FT8 Subtract)
    """

    def __init__(self):
        gr.basic_block.__init__(self, name='FT8 Demodulator', in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern('cand'))
        self.message_port_register_out(pmt.intern('llr'))
        self.message_port_register_out(pmt.intern('audio'))
        self.set_msg_handler(pmt.intern('cand'), self.handle)

    def handle(self, m):
        d = from_msg(m)
        x = np.asarray(d['audio'], dtype=np.float64)
        self.message_port_pub(pmt.intern('audio'), to_msg(dict(slot=d['slot'], **{'pass': d['pass']}, audio=d['audio'])))
        f1s, ibs, llrs, css = [], [], [], []
        if len(d['f0']):
            ds = Downsampler(x)
            for f0, xdt in zip(d['f0'], d['xdt']):
                f1, ib, cdp = fine_sync(ds, float(f0), float(xdt))
                cs = symbol_spectra(cdp, ib)
                hard = np.argmax(np.abs(cs), axis=1)
                if (hard[SYNC_SYMS] == _SYNC_TONE).sum() <= 6:
                    continue
                f1s.append(f1)
                ibs.append(ib)
                llrs.append(soft_metrics(cs))
                css.append(cs)
        out = {k: d[k] for k in ('slot', 'label', 'pretty', 'pass', 't0', 'base')}
        out['f1'] = np.array(f1s, dtype=np.float64)
        out['ib'] = np.array(ibs, dtype=np.int32)
        out['llr'] = np.array(llrs, dtype=np.float32).reshape(-1)
        out['cs'] = np.array(css, dtype=np.complex64).reshape(-1)
        self.message_port_pub(pmt.intern('llr'), to_msg(out))
