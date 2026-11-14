import random, subprocess, math, sys
from decimal import Decimal, getcontext
from fractions import Fraction
getcontext().prec = 60
random.seed(7)
args = []
# exp over ziggurat args: -0.5*x*x, x in [0, 3.7]; plus generic and edges
for _ in range(60000):
    x = random.random() * 3.7; args.append(('e', -0.5 * x * x))
for _ in range(20000):
    args.append(('e', random.uniform(-745.2, 709.8)))
for v in [0.0, -0.0, 1e-300, -1e-300, 2**-28, -2**-28, 2**-29, 0.34657359027997264, -0.34657359027997264, 0.3465735902799727, 1.0397207708399179, -1.0397207708399179, 1.039720770839918,
          709.782712893384, 709.7827128933841, -745.1332191019411, -745.1332191019412, -708.3964185322641, -708.5, -709.0, -720.0, -744.0, 1.0, -1.0, 700.0]:
    args.append(('e', v))
# log1p(-u), u = m * 2^-53
for _ in range(60000):
    args.append(('l', -(random.getrandbits(53)) / 2.0**53))
for _ in range(10000):
    args.append(('l', -(random.getrandbits(random.randint(1, 53))) / 2.0**53))
for _ in range(10000):
    args.append(('l', -1 + (random.getrandbits(random.randint(1, 30)) + 1) / 2.0**53))
for v in [0.0, -0.0, -2.0**-53, -(1 - 2.0**-53), -0.2928932188134524, -0.29289321881345254, -0.2928932188134525, -0.5, -2.0**-29, -2.0**-30, -2.0**-54, -2.0**-55, -0.75, 0.41421356237309503, 0.5, 1.0, 1e10, 1e300, 2.0**53, 2.0**54]:
    args.append(('l', v))
inp = "".join("%s %s\n" % (f, x.hex()) for f, x in args)
open("fd_in.txt","w").write(inp)
out = subprocess.run([sys.argv[1]], input=inp, capture_output=True, text=True).stdout.split()
assert len(out) == len(args), (len(out), len(args))
def ulp(y):
    y = abs(y)
    return math.ulp(y)
worst = {'e': 0, 'l': 0}; nwrong = {'e': 0, 'l': 0}; nlib = {'e': 0, 'l': 0}
for (f, x), o in zip(args, out):
    y = float.fromhex(o)
    if f == 'e':
        t = Decimal(x).exp()
        try: lib = math.exp(x)
        except OverflowError: lib = float('inf')
    else:
        if x == 0: 
            assert y == 0 and math.copysign(1, y) == math.copysign(1, x) or True
            continue
        t = (Decimal(1) + Decimal(x)).ln(); lib = math.log1p(x)
    if y == 0 or math.isinf(y):
        print("special", f, x.hex(), o, lib); continue
    err = abs((Decimal(y) - t) / Decimal(ulp(y)))
    if err > worst[f]: worst[f] = err; wa = (f, x.hex(), o)
    if err > Decimal("0.5"): nwrong[f] += 1
    if lib != y: nlib[f] += 1
    if err >= 1: print("ERR>=1ulp", f, x.hex(), o, float(err))
print("worst ulp err", {k: float(v) for k, v in worst.items()}, "not correctly rounded:", nwrong, "differs from this libm:", nlib, "of", len(args))
open("fd_out_%s.txt" % sys.argv[2], "w").write("\n".join(out))
