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
_TSTEP = 480                    # spectrogram hop: 1/4 symbol (40 ms)
_NFFT1 = 3840                   # 1920 samples zero-padded -> 3.125 Hz bins
_JSTRT = 12                     # 0.5 s / 40 ms
_JZ = 62                        # lag search +-62 steps (+-2.48 s)
_MLAG = 10                      # "narrow" lag search +-0.4 s


def sync_search(x, fmin, fmax, syncmin=1.3, maxcand=300):
    """Coarse Costas search. Returns list of (f0 Hz, xdt s, sync) and the baseline
    (40th-percentile) spectrum power per 3.125 Hz bin."""
    nhsym = len(x) // _TSTEP - 3
    idx = np.arange(nhsym)[:, None] * _TSTEP + np.arange(NSPS)[None, :]
    X = np.fft.rfft(x[idx], _NFFT1, axis=1)
    s = (X.real ** 2 + X.imag ** 2).T                     # (1921 bins, nhsym)
    df = FS / _NFFT1
    ia = max(1, int(round(fmin / df)))
    ib = min(int(round(fmax / df)), s.shape[0] - 15)
    nf = ib - ia + 1
    pad = np.zeros((s.shape[0], _JZ + nhsym + _JZ + 300))
    pad[:, _JZ:_JZ + nhsym] = s
    tones = np.stack([pad[ia + 2 * k: ia + 2 * k + nf] for k in range(7)])   # (7,nf,T)
    t_g = np.zeros((3, nf, 2 * _JZ + 1))
    t0_g = np.zeros((3, nf, 2 * _JZ + 1))
    for g in range(3):
        for n in range(7):
            m0 = _JSTRT + 4 * n + 144 * g                  # lag j=0 position
            sl = tones[:, :, m0: m0 + 2 * _JZ + 1]        # j = -JZ..JZ (pad offset JZ cancels)
            t_g[g] += sl[COSTAS[n]]
            t0_g[g] += sl.sum(axis=0)
    eps = 1e-30
    t = t_g.sum(0)
    t0 = (t0_g.sum(0) - t) / 6.0
    sync_abc = t / (t0 + eps)
    tbc = t_g[1] + t_g[2]
    t0bc = (t0_g[1] + t0_g[2] - tbc) / 6.0
    sync_bc = tbc / (t0bc + eps)
    sync2d = np.maximum(sync_abc, sync_bc)                 # (nf, 2JZ+1)
    narrow = sync2d[:, _JZ - _MLAG: _JZ + _MLAG + 1]
    jp1 = np.argmax(narrow, axis=1) - _MLAG
    red1 = narrow.max(axis=1)
    jp2 = np.argmax(sync2d, axis=1) - _JZ
    red2 = sync2d.max(axis=1)
    base1 = np.sort(red1)[int(round(0.40 * nf))]
    base2 = np.sort(red2)[int(round(0.40 * nf))]
    red1 = red1 / max(base1, eps)
    red2 = red2 / max(base2, eps)
    cand = []
    for red, jp in ((red1, jp1), (red2, jp2)):
        order = np.argsort(-red)[:maxcand]
        for i in order:
            if red[i] < syncmin:
                break
            cand.append(((ia + i) * df, (jp[i] - 0.5) * (_TSTEP / FS), float(red[i])))
    # remove near-duplicates (keep stronger)
    cand.sort(key=lambda c: -c[2])
    keep = []
    for c in cand:
        if all(abs(c[0] - k[0]) >= 4.0 or abs(c[1] - k[1]) >= 0.04 for k in keep):
            keep.append(c)
        if len(keep) >= maxcand:
            break
    # baseline noise spectrum (per 3.125 Hz bin), robust to signals
    base_spec = np.percentile(s, 40, axis=1)
    return keep, base_spec




# ================================================================ Stage 2 block
class blk(gr.basic_block):
    """FT8 Sync Search: slot audio -> candidate list.

    Builds a 3.125 Hz x 40 ms power spectrogram of the slot and scores the three
    7x7 Costas sync arrays at every frequency and time lag (-2.5 .. +2.5 s).
    Pass 0 audio comes from the Slot Clock; later passes are residual audio from
    FT8 Subtract (same input port).
    Output 'cand': input fields + f0[], xdt[], sync[], base[] (noise baseline per bin)
    """

    def __init__(self, f_min=200.0, f_max=3000.0, syncmin=1.2, maxcand=400):
        gr.basic_block.__init__(self, name='FT8 Sync Search', in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern('slot'))
        self.message_port_register_out(pmt.intern('cand'))
        self.set_msg_handler(pmt.intern('slot'), self.handle)
        self.f_min, self.f_max = float(f_min), float(f_max)
        self.syncmin, self.maxcand = float(syncmin), int(maxcand)

    def handle(self, m):
        d = from_msg(m)
        x = np.asarray(d['audio'], dtype=np.float64)
        if d['pass'] == 0:
            x = x / (np.std(x) + 1e-12)               # unit power; residuals keep this scale
        cands, base = sync_search(x, self.f_min, self.f_max, self.syncmin, self.maxcand) if np.any(x) else ([], np.zeros(1))
        d['audio'] = x.astype(np.float32)
        d['f0'] = np.array([c[0] for c in cands], dtype=np.float64)
        d['xdt'] = np.array([c[1] for c in cands], dtype=np.float64)
        d['sync'] = np.array([c[2] for c in cands], dtype=np.float64)
        if 'base' not in d:                            # keep the pass-0 noise baseline for SNR
            d['base'] = base.astype(np.float64)
        self.message_port_pub(pmt.intern('cand'), to_msg(d))
