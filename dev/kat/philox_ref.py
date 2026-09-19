"""Independent philox4x64-10, written from the Salmon et al. (2011) definition.
Validated against Random123's published kat_vectors before it is trusted to
generate anything."""
M0 = 0xD2E7470EE14C6C93   # PHILOX_M4x64_0
M1 = 0xCA5A826395121157   # PHILOX_M4x64_1
W0 = 0x9E3779B97F4A7C15   # golden ratio
W1 = 0xBB67AE8584CAA73B   # sqrt(3)-1
MASK = (1 << 64) - 1

def _round(c, k):
    p0 = M0 * c[0]; hi0, lo0 = p0 >> 64, p0 & MASK
    p1 = M1 * c[2]; hi1, lo1 = p1 >> 64, p1 & MASK
    return [hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0]

def philox4x64(ctr, key, rounds=10):
    c, k = list(ctr), list(key)
    for r in range(rounds):
        if r:                                  # bump before every round but the first
            k = [(k[0] + W0) & MASK, (k[1] + W1) & MASK]
        c = _round(c, k)
    return c

def hx(v): return " ".join("%016x" % x for x in v)

# ---- validate against upstream's published vectors -------------------------
KAT = [
 (10, [0]*4, [0]*2,
  "16554d9eca36314c db20fe9d672d0fdc d7e772cee186176b 7e68b68aec7ba23b"),
 (10, [MASK]*4, [MASK]*2,
  "87b092c3013fe90b 438c3c67be8d0224 9cc7d7c69cd777b6 a09caebf594f0ba0"),
 (10, [0x243f6a8885a308d3, 0x13198a2e03707344, 0xa4093822299f31d0, 0x082efa98ec4e6c89],
      [0x452821e638d01377, 0xbe5466cf34e90c6c],
  "a528f45403e61d95 38c72dbd566e9788 a5a1610e72fd18b5 57bd43b5e52b7fe6"),
 (7, [0]*4, [0]*2,
  "5dc8ee6268ec62cd 139bc570b6c125a0 84d6deb4fb65f49e aff7583376d378c2"),
 (7, [MASK]*4, [MASK]*2,
  "071dd84367903154 48e2bbdc722b37d1 6afa9890bb89f76c 9194c8d8ada56ac7"),
 (7, [0x243f6a8885a308d3, 0x13198a2e03707344, 0xa4093822299f31d0, 0x082efa98ec4e6c89],
      [0x452821e638d01377, 0xbe5466cf34e90c6c],
  "513a366704edf755 f05d9924c07044d3 bef2cb9cbea74c6c 8db948de4caa1f8a"),
]
ok = True
for rounds, ctr, key, want in KAT:
    got = hx(philox4x64(ctr, key, rounds))
    good = got == want
    ok &= good
    print(f"  rounds={rounds:2d}  {'OK  ' if good else 'FAIL'}  {got}")
print("\nreference implementation validated:", ok)
