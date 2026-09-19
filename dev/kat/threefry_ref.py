"""Independent threefry4x64, written from the Threefish/Skein round structure
Random123 implements. Validated against Random123's published kat_vectors
before it is trusted to generate anything."""
MASK = (1 << 64) - 1
KS_PARITY = 0x1BD11BDAA9FC1A22
ROT = [(14,16),(52,57),(23,40),(5,37),(25,33),(46,12),(58,22),(32,32)]

def _rotl(x, n): return ((x << n) | (x >> (64 - n))) & MASK

def threefry4x64(ctr, key, rounds=20):
    ks = list(key) + [KS_PARITY]
    for k in key:
        ks[4] ^= k
    X = [(ctr[i] + ks[i]) & MASK for i in range(4)]
    for r in range(rounds):
        r0, r1 = ROT[r % 8]
        if r % 2 == 0:                       # pair (0,1) and (2,3)
            X[0] = (X[0] + X[1]) & MASK; X[1] = _rotl(X[1], r0); X[1] ^= X[0]
            X[2] = (X[2] + X[3]) & MASK; X[3] = _rotl(X[3], r1); X[3] ^= X[2]
        else:                                # pair (0,3) and (2,1)
            X[0] = (X[0] + X[3]) & MASK; X[3] = _rotl(X[3], r0); X[3] ^= X[0]
            X[2] = (X[2] + X[1]) & MASK; X[1] = _rotl(X[1], r1); X[1] ^= X[2]
        if r % 4 == 3:                       # inject subkey every 4 rounds
            i = (r + 1) // 4
            for j in range(4):
                X[j] = (X[j] + ks[(i + j) % 5]) & MASK
            X[3] = (X[3] + i) & MASK
    return X

def hx(v): return " ".join("%016x" % x for x in v)

# Self-test. Guarded, because importing this module must not run it -- an
# unguarded sys.exit() here silently terminated emit_kat.py mid-run.
if __name__ == "__main__":
    PI = [0x243f6a8885a308d3, 0x13198a2e03707344, 0xa4093822299f31d0, 0x082efa98ec4e6c89]
    PK = [0x452821e638d01377, 0xbe5466cf34e90c6c, 0xc0ac29b7c97c50dd, 0x3f84d5b5b5470917]
    KAT = [
     (13,[0]*4,[0]*4,"4071fabee1dc8e05 02ed3113695c9c62 397311b5b89f9d49 e21292c3258024bc"),
     (13,[MASK]*4,[MASK]*4,"7eaed935479722b5 90994358c429f31c 496381083e07a75b 627ed0d746821121"),
     (13,PI,PK,"4361288ef9c1900c 8717291521782833 0d19db18c20cf47e a0b41d63ac8581e5"),
     (20,[0]*4,[0]*4,None), (72,[0]*4,[0]*4,None), (72,PI,PK,None),
    ]
    import sys, re
    # pull the 20- and 72-round expectations straight from the upstream file
    kv = open("/Users/pbtz/Documents/repos/gh/public/random123/tests/kat_vectors").read().split("\n")
    def upstream(rounds, ctr, key):
        for ln in kv:
            f = ln.split()
            if len(f) == 14 and f[0] == "threefry4x64" and int(f[1]) == rounds:
                if [int(x,16) for x in f[2:6]] == ctr and [int(x,16) for x in f[6:10]] == key:
                    return " ".join(f[10:14])
        return None
    ok = True
    for rounds, ctr, key, want in KAT:
        want = want or upstream(rounds, ctr, key)
        if want is None:                      # not published for this combination
            print(f"  rounds={rounds:2d}  skip  (no upstream vector)"); continue
        got = hx(threefry4x64(ctr, key, rounds))
        good = (got == want); ok &= good
        print(f"  rounds={rounds:2d}  {'OK  ' if good else 'FAIL'}  {got}")
        if not good: print(f"            want  {want}")
    print("\nreference validated:", ok)
    sys.exit(0 if ok else 1)
