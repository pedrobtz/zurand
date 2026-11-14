import sys, struct, array
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2]]
SRC, D = sys.argv[1], sys.argv[2]
exec(open(__file__.replace("compare.py", "oracle.py")).read())
T = load_tables(SRC)
key = splitmix_keys(42, 3)[1]
print("key", "%016x %016x" % key)
def rd(fn, code):
    a = array.array(code); b = open(fn, "rb").read(); a.frombytes(b); return a
def bits(x): return struct.unpack("<Q", struct.pack("<d", x))[0]
N = 300000
ranges = [(1,6), (-1, 2147483647), (-2147483647, 2147483647), (5,5), (-1610612736, 1610612736), (0, 2**31-2)]
offs = [511, 4095, 2**32 - 7, 2**53 - 20000]
bad = 0
def cmp(name, got, want, isf=True):
    global bad
    mm = [i for i in range(len(want)) if (bits(got[i]) != bits(want[i]) if isf else got[i] != want[i])]
    print("%-28s n=%d mismatches=%d" % (name, len(want), len(mm)), mm[:5], [(got[i], want[i]) for i in mm[:3]])
    bad += len(mm)
for e in ["philox4x64", "threefry4x64", "xoshiro256pp"]:
    f = lambda nm: D + "/" + e + "_" + nm
    hb = open(f("bits")).read().split()
    wb = words(e, key, 0, 0, len(hb))
    cmp(e + " bits", [int(h, 16) for h in hb], wb, False)
    w = words(e, key, 1, 0, N); cmp(e + " unif", rd(f("unif"), "d"), [u01(x) for x in w])
    w = words(e, key, 2, 0, N); cmp(e + " norm", rd(f("norm"), "d"), [zig(T, e, key, i, x) for i, x in enumerate(w)])
    w = words(e, key, 7, 0, N); cmp(e + " mcf", rd(f("mcf"), "d"), [mcf(T, e, key, i, x) for i, x in enumerate(w)])
    w = words(e, key, 3, 0, N)
    for r, (lo, hi) in enumerate(ranges):
        cmp(e + " int%d" % (r + 1), rd(f("int%d" % (r + 1)), "i"), [lemire(e, key, i, x, lo, hi) for i, x in enumerate(w)], False)
    for o, off in enumerate(offs):
        n = 20000
        w = words(e, key, 1, off, n); cmp(e + " unif@%d" % off, rd(f("unifoff%d" % (o + 1)), "d"), [u01(x) for x in w])
        w = words(e, key, 2, off, n); cmp(e + " norm@%d" % off, rd(f("normoff%d" % (o + 1)), "d"), [zig(T, e, key, off + i, x) for i, x in enumerate(w)])
        w = words(e, key, 7, off, n); cmp(e + " mcf@%d" % off, rd(f("mcfoff%d" % (o + 1)), "d"), [mcf(T, e, key, off + i, x) for i, x in enumerate(w)])
    print(stats)
print("TOTAL BAD", bad)
