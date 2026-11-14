"""From-spec oracle for zurand transforms (hunt-3). Pure Python, no package code."""
import re, struct, math, sys
MASK = (1 << 64) - 1
SRC = sys.argv[1] if len(sys.argv) > 1 else None
M0 = 0xD2E7470EE14C6C93; M1 = 0xCA5A826395121157
W0 = 0x9E3779B97F4A7C15; W1 = 0xBB67AE8584CAA73B
def philox(ctr, key, rounds=10):
    c0, c1, c2, c3 = ctr; k0, k1 = key
    for r in range(rounds):
        if r:
            k0 = (k0 + W0) & MASK; k1 = (k1 + W1) & MASK
        p0 = M0 * c0; p1 = M1 * c2
        c0, c1, c2, c3 = (p1 >> 64) ^ c1 ^ k0, p1 & MASK, (p0 >> 64) ^ c3 ^ k1, p0 & MASK
    return [c0, c1, c2, c3]
assert " ".join("%016x" % x for x in philox([0]*4, [0]*2)) == "16554d9eca36314c db20fe9d672d0fdc d7e772cee186176b 7e68b68aec7ba23b"
KS_PARITY = 0x1BD11BDAA9FC1A22
ROT = [(14,16),(52,57),(23,40),(5,37),(25,33),(46,12),(58,22),(32,32)]
def _rotl(x, n): return ((x << n) | (x >> (64 - n))) & MASK
def threefry(ctr, key, rounds=13):
    ks = list(key) + [KS_PARITY]
    for k in key: ks[4] ^= k
    X = [(ctr[i] + ks[i]) & MASK for i in range(4)]
    for r in range(rounds):
        r0, r1 = ROT[r % 8]
        if r % 2 == 0:
            X[0] = (X[0] + X[1]) & MASK; X[1] = _rotl(X[1], r0); X[1] ^= X[0]
            X[2] = (X[2] + X[3]) & MASK; X[3] = _rotl(X[3], r1); X[3] ^= X[2]
        else:
            X[0] = (X[0] + X[3]) & MASK; X[3] = _rotl(X[3], r0); X[3] ^= X[0]
            X[2] = (X[2] + X[1]) & MASK; X[1] = _rotl(X[1], r1); X[1] ^= X[2]
        if r % 4 == 3:
            i = (r + 1) // 4
            for j in range(4): X[j] = (X[j] + ks[(i + j) % 5]) & MASK
            X[3] = (X[3] + i) & MASK
    return X
assert " ".join("%016x" % x for x in threefry([0]*4, [0]*4)) == "4071fabee1dc8e05 02ed3113695c9c62 397311b5b89f9d49 e21292c3258024bc"

def splitmix_keys(seed, n):
    s = seed; out = []
    def nxt():
        nonlocal s
        s = (s + 0x9e3779b97f4a7c15) & MASK
        z = s
        z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & MASK
        z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & MASK
        return z ^ (z >> 31)
    for _ in range(n):
        a = nxt(); b = nxt(); out.append((a, b))
    return out

def block(engine, key, c0, c1, purpose):
    if engine == "threefry4x64":
        k = [key[0], key[1], key[0] ^ 0x9e3779b97f4a7c15, key[1] ^ 0xbb67ae8584caa73b]
        return threefry([c0, c1, purpose, 0], k)
    return philox([c0, c1, purpose, 1 if engine == "xoshiro256pp" else 0], key)

def words(engine, key, purpose, start, n):
    """words at positions start .. start+n-1"""
    out = []
    if engine == "xoshiro256pp":
        pos = start
        while len(out) < n:
            sub, off = divmod(pos, 512)
            st = block(engine, key, sub, 0, purpose)
            if not any(st): st[0] = 1
            s0, s1, s2, s3 = st
            take = min(512 - off, n - len(out))
            for j in range(off + take):
                r = (_rotl((s0 + s3) & MASK, 23) + s0) & MASK
                t = (s1 << 17) & MASK
                s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t
                s3 = _rotl(s3, 45)
                if j >= off: out.append(r)
            pos += take
    else:
        b = start >> 2
        blk = None
        for p in range(start, start + n):
            if blk is None or (p >> 2) != b:
                b = p >> 2; blk = block(engine, key, b, 0, purpose)
            out.append(blk[p & 3])
    return out

class Retry:
    def __init__(self, engine, key, index, purpose):
        self.e, self.k, self.i, self.p = engine, key, index, purpose
        self.g = 0; self.buf = []
    def next(self):
        if not self.buf:
            self.g += 1
            self.buf = block(self.e, self.k, self.i, self.g, self.p)
        return self.buf.pop(0)

def u01(w):   # (m + 0.5) * 2^-52
    return ((w >> 12) * 2 + 1) / 9007199254740992.0
def u53(w): return (w >> 11) * (1.0 / 9007199254740992.0)

# ---- tables ----
def load_tables(src):
    t = open(src + "/src/numpyzig/ziggurat_constants.h").read()
    def arr(name):
        m = re.search(name + r"\[\]\s*=\s*\{(.*?)\};", t, re.S); return m.group(1)
    ki = [int(x[:-3], 16) for x in re.findall(r"0x[0-9A-Fa-f]+ULL", arr("ki_double"))]
    wi = [float(x) for x in re.findall(r"[-+0-9.eE]+", arr("wi_double")) if x.strip()]
    fi = [float(x) for x in re.findall(r"[-+0-9.eE]+", arr("fi_double")) if x.strip()]
    assert len(ki) == 256 and len(wi) == 256 and len(fi) == 256, (len(ki), len(wi), len(fi))
    r = float(re.search(r"ziggurat_nor_r = ([0-9.]+)", t).group(1))
    inv_r = float(re.search(r"ziggurat_nor_inv_r =\s*([0-9.]+)", t).group(1))
    m = open(src + "/src/mcfarland_normal.h").read()
    def harr(name):
        return re.search(name + r"\[\d+\]\s*=\s*\{(.*?)\};", m, re.S).group(1)
    X = [float.fromhex(x) for x in re.findall(r"0x[0-9a-fA-F.]+p[-+]\d+", harr("zurand_mcf_x"))]
    Y = [float.fromhex(x) for x in re.findall(r"0x[0-9a-fA-F.]+p[-+]\d+", harr("zurand_mcf_y"))]
    AT = [int(x, 16) for x in re.findall(r"UINT64_C\((0x[0-9a-f]+)\)", harr("zurand_mcf_alias_t"))]
    AJ = [int(x) for x in re.findall(r"\d+", harr("zurand_mcf_alias_j"))]
    x0 = float.fromhex(re.search(r"zurand_mcf_x0 = (\S+);", m).group(1))
    ix0 = float.fromhex(re.search(r"zurand_mcf_inv_x0 = (\S+);", m).group(1))
    assert len(X) == 254 and len(Y) == 254 and len(AT) == 256 and len(AJ) == 256
    return dict(ki=ki, wi=wi, fi=fi, r=r, inv_r=inv_r, X=X, Y=Y, AT=AT, AJ=AJ, x0=x0, ix0=ix0)

stats = dict(wedge=0, tail=0, tail_retry=0, wedge_retry=0, mcf_edge=0, mcf_tail=0, mcf_retry=0, lem_retry=0)
EXP = math.exp; LOG1P = math.log1p
def zig(T, engine, key, index, w):
    ki, wi, fi = T["ki"], T["wi"], T["fi"]
    idx = w & 0xff; sign = (w >> 8) & 1; rabs = (w >> 9) & 0x000fffffffffffff
    x = rabs * wi[idx]
    if sign: x = -x
    if rabs < ki[idx]:
        return x if rabs else 0.0
    s = Retry(engine, key, index, 2)
    first = True
    while True:
        Yw = s.next()
        if idx == 0:
            stats["tail"] += 1
            yy = -LOG1P(-u53(Yw))
            while True:
                xx = -T["inv_r"] * LOG1P(-u53(s.next()))
                if yy + yy > xx * xx:
                    return -(T["r"] + xx) if sign else T["r"] + xx
                stats["tail_retry"] += 1
                yy = -LOG1P(-u53(s.next()))
        stats["wedge"] += 1
        u = u53(Yw)
        if (fi[idx - 1] - fi[idx]) * u + fi[idx] < EXP(-0.5 * x * x):
            return x
        stats["wedge_retry"] += 1
        w = s.next()
        idx = w & 0xff; sign = (w >> 8) & 1; rabs = (w >> 9) & 0x000fffffffffffff
        x = rabs * wi[idx]
        if sign: x = -x
        if rabs < ki[idx]:
            return x

def s64(w): return w - (1 << 64) if w >> 63 else w
def mcf(T, engine, key, index, w):
    X, Y = T["X"], T["Y"]
    i = w & 0xff
    if i < 253:
        return X[i] * float(s64(w))
    stats["mcf_edge"] += 1
    s = Retry(engine, key, index, 7)
    sign = -1.0 if w >> 63 else 1.0
    a = s.next()
    j = a & 0xff
    if (a >> 8) >= T["AT"][j]: j = T["AJ"][j]
    def interp(v, j, u): return v[j] * 9223372036854775808.0 + float(u) * (v[j - 1] - v[j])
    if j == 0:
        stats["mcf_tail"] += 1
        while True:
            xx = -T["ix0"] * LOG1P(-u53(s.next()))
            yy = -LOG1P(-u53(s.next()))
            if yy + yy > xx * xx:
                return sign * (T["x0"] + xx)
            stats["mcf_retry"] += 1
    u1 = w & 0x7fffffffffffffff
    while True:
        u2 = s.next() >> 1
        if j < 204 and u2 < u1: u1, u2 = u2, u1
        x = interp(X, j, u1)
        if interp(Y, j, u2) < EXP(-0.5 * x * x):
            return sign * x
        stats["mcf_retry"] += 1
        u1 = s.next() >> 1

def lemire(engine, key, index, w, lo, hi):
    rng = hi - lo + 1
    thr = (1 << 32) % rng
    x = w & 0xffffffff
    g = 0
    while True:
        p = x * rng
        if (p & 0xffffffff) >= thr:
            return lo + (p >> 32)
        stats["lem_retry"] += 1
        g += 1
        x = block(engine, key, index >> 2, g, 3)[index & 3] & 0xffffffff
