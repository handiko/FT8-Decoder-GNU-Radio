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


# ---------------------------------------------------------------- signal subtraction
def _gfsk_pulse(bt, t):
    from math import erf, pi, sqrt, log
    c = pi * sqrt(2.0 / log(2.0))
    v = np.vectorize(erf)
    return 0.5 * (v(c * bt * (t + 0.5)) - v(c * bt * (t - 0.5)))


_PULSE = _gfsk_pulse(2.0, (np.arange(3 * NSPS) - 1.5 * NSPS) / NSPS)


def gfsk_ref(tones, f0):
    """Unit-amplitude complex FT8 waveform (79 symbols, GFSK BT=2) at f0 Hz, 12 kHz."""
    hmod = 2.0 * np.pi / NSPS
    dphi = np.zeros((NN + 2) * NSPS)
    for j, tn in enumerate(tones):
        dphi[j * NSPS: j * NSPS + 3 * NSPS] += hmod * tn * _PULSE
    dphi[0:2 * NSPS] += hmod * tones[0] * _PULSE[NSPS:3 * NSPS]
    dphi[NN * NSPS:(NN + 2) * NSPS] += hmod * tones[-1] * _PULSE[0:2 * NSPS]
    dphi += 2.0 * np.pi * f0 / FS
    phi = np.concatenate([[0.0], np.cumsum(dphi[NSPS:NSPS + NN * NSPS - 1])])
    return np.exp(1j * phi)


_NFILT = 4000
_WFILT = np.cos(np.pi * (np.arange(_NFILT) - _NFILT // 2) / _NFILT) ** 2
_WFILT /= _WFILT.sum()


_DEN_CACHE = {}


def _fftconv_same(v):
    n = len(v) + _NFILT - 1
    nf = 1 << int(np.ceil(np.log2(n)))
    r = np.fft.ifft(np.fft.fft(v, nf) * np.fft.fft(_WFILT, nf))[:n]
    s = (_NFILT - 1) // 2
    return r[s:s + len(v)]


def subtract(x, tones, f0, t0):
    """Remove one decoded signal from x (in place). t0 = start time in s from slot start."""
    cref = gfsk_ref(tones, f0)
    n0 = int(round(t0 * FS))
    a, b = max(0, n0), min(len(x), n0 + len(cref))
    if b - a < NSPS:
        return
    cr = cref[a - n0: b - n0]
    camp = x[a:b] * np.conj(cr)
    num = _fftconv_same(camp)
    den = _DEN_CACHE.get(b - a)
    if den is None:                                   # edge correction, same for equal lengths
        den = np.maximum(_fftconv_same(np.ones(b - a)).real, 1e-6)
        _DEN_CACHE[b - a] = den
    cfilt = num / den
    x[a:b] -= 2.0 * np.real(cfilt * cr)


# ================================================================ Stage 6 block
class blk(gr.basic_block):
    """FT8 Subtract: remove decoded signals and send the residual back to the search.

    Each new decode is rebuilt as an ideal FT8 waveform (GFSK, BT = 2) at its
    frequency and start time; its slowly varying amplitude and phase are measured
    with a ~0.33 s smoothing filter, and the fitted signal is subtracted. The residual
    goes back to FT8 Sync Search as the next pass, where weaker signals that were
    hidden under the strong ones can now be found.
    In 'audio': the audio each pass was decoded from (from FT8 Demodulator)
    In 'new'  : the decodes that pass produced (from FT8 Unpack + Report)
    Out 'slot': residual audio, pass + 1
    """

    def __init__(self):
        gr.basic_block.__init__(self, name='FT8 Subtract', in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern('audio'))
        self.message_port_register_in(pmt.intern('new'))
        self.message_port_register_out(pmt.intern('slot'))
        self.set_msg_handler(pmt.intern('audio'), self.on_audio)
        self.set_msg_handler(pmt.intern('new'), self.on_new)
        self.cache = {}

    def on_audio(self, m):
        d = from_msg(m)
        self.cache[(d['slot'], d['pass'])] = np.asarray(d['audio'], dtype=np.float64)
        if len(self.cache) > 12:                      # forget audio of finished slots
            for k in sorted(self.cache)[:len(self.cache) - 12]:
                del self.cache[k]

    def on_new(self, m):
        d = from_msg(m)
        x = self.cache.pop((d['slot'], d['pass']), None)
        if x is None:
            return
        x = x.copy()
        tones = np.asarray(d['tones'], dtype=int).reshape(-1, NN)
        for t, f1, t0 in zip(tones, d['f1'], d['tstart']):
            subtract(x, t, float(f1), float(t0))
        self.message_port_pub(pmt.intern('slot'), to_msg(dict(
            slot=d['slot'], label=d['label'], pretty=d['pretty'], **{'pass': d['pass'] + 1}, t0=d['t0'],
            base=np.asarray(d['base'], dtype=np.float64), audio=x.astype(np.float32))))
