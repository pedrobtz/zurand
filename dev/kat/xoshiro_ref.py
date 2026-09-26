"""Independent xoshiro256++, from Vigna's published definition. Validated
against the author's own reference C (prng.di.unimi.it/xoshiro256plusplus.c,
public domain) run from a fixed state, before it is trusted here."""
MASK = (1 << 64) - 1
def _rotl(x, k): return ((x << k) | (x >> (64 - k))) & MASK

def xoshiro256pp(state, n):
    s0, s1, s2, s3 = state
    out = []
    for _ in range(n):
        r = (_rotl((s0 + s3) & MASK, 23) + s0) & MASK
        t = (s1 << 17) & MASK
        s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t
        s3 = _rotl(s3, 45)
        out.append(r)
    return out

# zurand's chunk seeding: Philox at {chunk, 0, purpose, 1} gives the state,
# with the all-zero guard. The final 1 is the xoshiro engine tag in counter
# word 3, which keeps these seeds apart from the philox engine's own blocks.
XOSHIRO_TAG = 1
def zurand_xoshiro_chunk(philox4x64, key, chunk, purpose, nwords):
    st = philox4x64([chunk, 0, purpose, XOSHIRO_TAG], key, 10)
    if (st[0] | st[1] | st[2] | st[3]) == 0: st[0] = 1
    return xoshiro256pp(st, nwords)

if __name__ == "__main__":
    import sys
    ref = [int(l, 16) for l in open("/tmp/xo_ref_out.txt").read().split()]
    got = xoshiro256pp([0x0123456789abcdef, 0xfedcba9876543210,
                        0xdeadbeefcafebabe, 0x0000000000000001], len(ref))
    ok = got == ref
    for g, r in zip(got, ref):
        print(f"  {'OK  ' if g == r else 'FAIL'}  {g:016x}")
    print("\nreference validated against Vigna's C:", ok)
    sys.exit(0 if ok else 1)
