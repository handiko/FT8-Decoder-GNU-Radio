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


# ---------------------------------------------------------------- LDPC(174,91) parity checks
_NM_TXT = (
    '4,31,59,91,92,96,153;5,32,60,93,115,146;6,24,61,94,122,151;7,33,62,95,96,143;8,25,63,83,93,96,148;6,32,64,97,126,138;5,34,65,78,98,107,154;9,35,66,99,139,146;'
    '10,36,67,100,107,126;11,37,67,87,101,139,158;12,38,68,102,105,155;13,39,69,103,149,162;8,40,70,82,104,114,145;14,41,71,88,102,123,156;15,42,59,106,123,159;1,33,72,106,107,157;'
    '16,43,73,108,141,160;17,37,74,81,109,131,154;11,44,75,110,121,166;45,55,64,111,130,161,173;8,46,71,112,119,166;18,36,76,89,113,114,143;19,38,77,104,116,163;20,47,70,92,138,165;'
    '2,48,74,113,128,160;21,45,78,83,117,121,151;22,47,58,118,127,164;16,39,62,112,134,158;23,43,79,120,131,145;19,35,59,73,110,125,161;20,36,63,94,136,161;14,31,79,98,132,164;'
    '3,44,80,124,127,169;19,46,81,117,135,167;7,49,58,90,100,105,168;12,50,61,118,119,144;13,51,64,114,118,157;24,52,76,129,148,149;25,53,69,90,101,130,156;20,46,65,80,120,140,170;'
    '21,54,77,100,140,171;35,82,133,142,171,174;14,30,83,113,125,170;4,29,68,120,134,173;1,4,52,57,86,136,152;26,51,56,91,122,137,168;52,84,110,115,145,168;7,50,81,99,132,173;'
    '23,55,67,95,172,174;26,41,77,109,141,148;2,27,41,61,62,115,133;27,40,56,124,125,126;18,49,55,124,141,167;6,33,85,108,116,156;28,48,70,85,105,129,158;9,54,63,131,147,155;'
    '22,53,68,109,121,174;3,13,48,78,95,123;31,69,133,150,155,169;12,43,66,89,97,135,159;5,39,75,102,136,167;2,54,86,101,135,164;15,56,87,108,119,171;10,44,82,91,111,144,149;'
    '23,34,71,94,127,153;11,49,88,92,142,157;29,34,87,97,147,162;30,50,60,86,137,142,162;10,53,66,84,112,128,165;22,57,85,93,140,159;28,32,72,103,132,166;28,29,84,88,117,143,150;'
    '1,26,45,80,128,147;17,27,89,103,116,153;51,57,98,163,165,172;21,37,73,138,152,169;16,47,76,130,137,154;3,24,30,72,104,139;9,40,90,106,134,151;15,58,60,74,111,150,163;'
    '18,42,79,144,146,152;25,38,65,99,122,160;17,42,75,129,170,172'
)

_NM = [[int(v) - 1 for v in r.split(',')] for r in _NM_TXT.split(';')]

# Tanner graph as edge lists
_EC, _EV = [], []
for _c, _row in enumerate(_NM):
    for _v in _row:
        _EC.append(_c)
        _EV.append(_v)
_EC = np.array(_EC)
_EV = np.array(_EV)
_NE = len(_EC)                                       # 522 edges
_CHK_E = np.full((83, 7), _NE)                       # edge ids per check, pad -> dummy edge
_VAR_E = np.zeros((174, 3), dtype=int)               # edge ids per variable (3 each)
_fill_c = np.zeros(83, dtype=int)
_fill_v = np.zeros(174, dtype=int)
for _e in range(_NE):
    _c, _v = _EC[_e], _EV[_e]
    _CHK_E[_c, _fill_c[_c]] = _e
    _fill_c[_c] += 1
    _VAR_E[_v, _fill_v[_v]] = _e
    _fill_v[_v] += 1
_CHK_V = np.where(_CHK_E < _NE, _EV[np.minimum(_CHK_E, _NE - 1)], 174)  # var per check slot, pad -> 174


# ---------------------------------------------------------------- FEC: CRC-14, LDPC
def crc14(bits):
    """CRC-14 (poly 0x2757) of rows of bits, shape (B, n) -> (B,) ints."""
    bits = np.atleast_2d(bits).astype(np.int64)
    reg = np.zeros(bits.shape[0], dtype=np.int64)
    for k in range(bits.shape[1]):
        reg ^= bits[:, k] << 13
        top = (reg >> 13) & 1
        reg = ((reg << 1) ^ (top * 0x2757)) & 0x3FFF
    return reg


def crc_ok(a91):
    """a91: (B, 91) bits = 77 payload + 14 CRC. Valid and not all-zero."""
    msg = np.concatenate([a91[:, :77], np.zeros((a91.shape[0], 5), dtype=a91.dtype)], axis=1)
    rx = np.zeros(a91.shape[0], dtype=np.int64)
    for k in range(14):
        rx = (rx << 1) | a91[:, 77 + k].astype(np.int64)
    return (crc14(msg) == rx) & (a91[:, :77].sum(axis=1) > 0)


def ldpc_bp(llr, max_iter=30):
    """Sum-product BP. llr: (B,174), positive = bit 1.
    Returns (bits (B,174) uint8, ok (B,) bool: all parity checks satisfied)."""
    B = llr.shape[0]
    L = -np.asarray(llr, dtype=np.float64)           # work with log(P0/P1)
    m_cv = np.zeros((B, _NE + 1))                    # +1: dummy edge
    bits = np.zeros((B, 174), dtype=np.uint8)
    done = np.zeros(B, dtype=bool)
    out = np.zeros((B, 174), dtype=np.uint8)
    for _ in range(max_iter):
        tot = L + m_cv[:, _VAR_E].sum(axis=2)        # (B,174)
        bits = (tot < 0).astype(np.uint8)
        syn = bits[:, np.minimum(_CHK_V, 173)] * (_CHK_V < 174)
        ok = (syn.sum(axis=2) % 2 == 0).all(axis=1)
        newly = ok & ~done
        out[newly] = bits[newly]
        done |= ok
        if done.all():
            break
        m_vc = tot[:, _EV] - m_cv[:, :_NE]           # (B,E)
        t = np.tanh(np.clip(m_vc, -40, 40) / 2.0)
        t = np.concatenate([t, np.ones((B, 1))], axis=1)
        tc = t[:, _CHK_E]                            # (B,83,7)
        sgn = np.where(tc < 0, -1.0, 1.0)
        la = np.log(np.maximum(np.abs(tc), 1e-30))
        prod_l = la.sum(axis=2, keepdims=True) - la
        prod_s = sgn.prod(axis=2, keepdims=True) * sgn
        v = np.clip(prod_s * np.exp(prod_l), -0.999999, 0.999999)
        msg = 2.0 * np.arctanh(v)
        m_new = np.zeros((B, _NE + 1))
        m_new[:, _CHK_E.ravel()] = msg.reshape(B, -1)
        m_new[:, _NE] = 0.0
        m_cv = m_new
    out[~done] = bits[~done]
    return out, done


# ================================================================ Stage 4 block
class blk(gr.basic_block):
    """FT8 LDPC + CRC: soft bits -> codewords.

    Sum-product belief propagation on the (174,91) LDPC code for every candidate and
    all four soft-bit metrics at once (one numpy batch). A candidate is accepted when
    the parity checks are met, the CRC-14 matches, the payload is not all-zero and at
    most 36 hard-decision bits were corrected; the first metric that passes is used.
    Output 'cw': slot info + cw[K*174], f1[], ib[], cs[K*79*8], base[]
    """

    def __init__(self, max_iter=40):
        gr.basic_block.__init__(self, name='FT8 LDPC + CRC', in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern('llr'))
        self.message_port_register_out(pmt.intern('cw'))
        self.set_msg_handler(pmt.intern('llr'), self.handle)
        self.max_iter = int(max_iter)

    def handle(self, m):
        d = from_msg(m)
        n = len(d['f1'])
        keep, cws = [], []
        if n:
            LL = np.asarray(d['llr'], dtype=np.float64).reshape(n * 4, 174)
            bits, ok = ldpc_bp(LL, self.max_iter)
            good = ok & crc_ok(bits[:, :91])
            good &= ((LL > 0).astype(np.uint8) != bits).sum(axis=1) <= 36
            for c in range(n):
                rows = [4 * c + p for p in range(4) if good[4 * c + p]]
                if rows:
                    keep.append(c)
                    cws.append(bits[rows[0]])
        cs = np.asarray(d['cs']).reshape(n, NN, 8) if n else np.zeros((0, NN, 8), np.complex64)
        out = {k: d[k] for k in ('slot', 'label', 'pretty', 'pass', 't0', 'base')}
        out['cw'] = np.array(cws, dtype=np.uint8).reshape(-1)
        out['f1'] = np.asarray(d['f1'])[keep].astype(np.float64)
        out['ib'] = np.asarray(d['ib'])[keep].astype(np.int32)
        out['cs'] = cs[keep].reshape(-1)
        self.message_port_pub(pmt.intern('cw'), to_msg(out))
