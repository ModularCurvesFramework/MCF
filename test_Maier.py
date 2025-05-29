import pytest

from sage.all import QQ
from parameters import *

# Table 6 of https://arxiv.org/pdf/math/0611041.pdf
# {(N, N/d): (t_d, t'_d)}
# (N, ell) maps X0(N) to X0(N/ell)
# => ! Need to swap d and ell
MAIER_COVERS = {
    # 16 = 2^4
    (2, 2): (lambda t: (t+16)**3/t, lambda t: (t+256)**3/t**2),
    (4, 2): (lambda t: t*(t+16), lambda t: t**2/(t+16)),
    (8, 2): (lambda t: t*(t+8), lambda t: t**2/(t+4)),
    (16, 2): (lambda t: t*(t+4), lambda t: t**2/(t+2)),

    # 6 = 2*3
    (6, 2): (lambda t: t*(t+9)**2/(t+8), lambda t: t**2*(t+9)/(t+8)**2),
    (6, 3): (lambda t: t*(t+8)**3/(t+9), lambda t: t**3*(t+8)/(t+9)**3),

    # 9 = 3^2
    (3, 3): (lambda t: (t+27)*(t+3)**3/t, lambda t: (t+27)*(t+243)**3/t**3),
    (9, 3): (lambda t: t*(t**2+9*t+27), lambda t: t**3/(t**2+9*t+27)),

    # 12 = 2^2 * 3
    (12, 2): (lambda t: t*(t+6), lambda t: t**2/(t+2)),
    (12, 3): (lambda t: t*(t+4)**3/(t+3), lambda t: t**3*(t+4)/(t+3)**3),

    # 18 = 2 * 3^2
    (18, 2): (lambda t: t*(t+3)**2/(t+2), lambda t: t**2*(t+3)/(t+2)**2),
    (18, 3): (lambda t: t*(t**2+6*t+12), lambda t: t**3/(t**2+3*t+3)),

    # 25 = 5^2
    (5, 5): (lambda t: (t**2+10*t+5)**3/t, lambda t: (t**2+250*t+3125)**3/t**5),
    (25, 5): (lambda t: t*(t**4+5*t**3+15*t**2+25*t+25), lambda t: t**5/(t**4+5*t**3+15*t**2+25*t+25)),

    # 10 = 2*5
    (10, 2): (lambda t: t*(t+5)**2/(t+4), lambda t: t**2*(t+5)/(t+4)**2),
    (10, 5): (lambda t: t*(t+4)**5/(t+5), lambda t: t**5*(t+4)/(t+5)**5),

    # 7
    (7, 7): (lambda t: (t**2+13*t+49)*(t**2+5*t+1)**3/t, lambda t: (t**2+13*t+49)*(t**2+245*t+2401)**3/t**7),

    # 13
    (13, 13): (lambda t: (t**2+5*t+13)*(t**4+7*t**3+20*t**2+19*t+1)**3/t, lambda t: (t**2+5*t+13)*(t**4+247*t**3+3380*t**2+15379*t+28561)**3/t**13),
}
MAIER_COVERS_TUPLES = [
    k + v for k, v in MAIER_COVERS.items()
]

def invert_linear(f):
    b = f(0)
    a = f(1) - b
    assert f(100) == a*100 + b
    assert f(101) == a*101 + b
    ret = lambda t: (t - b) / a
    assert ret(f(100)) == 100
    assert ret(f(123)) == 123
    return ret


# old test, includes derivation of inverse maps
PROJ = {
    2:  lambda t: t,
    4: lambda t: -4*t-8,
    8: lambda t: 4*t-4,
    16: lambda t: 2*t-2,
    6: lambda t: t-9,
    18: lambda t: t-2,
    12: lambda t: t-3,
    9: lambda t: t-3,
    25: lambda t: t-1,
    10: lambda t: t-4,
}
iPROJ = {
    N: invert_linear(proj) for N, proj in PROJ.items()
}

t = QQ['t'].gen()
TO_MAIER = {N: proj(t) for N, proj in PROJ.items()}
FROM_MAIER = {N: proj(t) for N, proj in iPROJ.items()}

for N in (1, 2, 3, 5, 7, 13):
    FROM_MAIER[N] = t
    TO_MAIER[N] = t


@pytest.mark.parametrize('N,l,Maier_left,Maier_right', MAIER_COVERS_TUPLES)
def test_Maier1(N, l, Maier_left, Maier_right):
    L = globals()["Level%d" % N]

    Maier_left = Maier_left(t)
    Maier_right = Maier_right(t)
    t_my = FROM_MAIER[N](t)

    # note: left/right swapped compared to Maier
    left = L(t_my).left(l).value
    left2 = TO_MAIER[N//l](left)
    assert left2 == Maier_right, (N, l, left2, Maier_right)

    right = L(t_my).right(l).value
    right2 = TO_MAIER[N//l](right)
    assert right2 == Maier_left, (N, l, right2, Maier_left)


@pytest.mark.parametrize('N,l,Maier_left,Maier_right', MAIER_COVERS_TUPLES)
def test_Maier2(N, l, Maier_left, Maier_right):
    # same test, different API (hardcoded maps in the parameters classes)
    t = QQ['t'].gen()
    L = globals()["Level%d" % N]

    Maier_left = Maier_left(t)
    Maier_right = Maier_right(t)

    # note: left/right swapped compared to Maier
    left2 = L.from_Maier(t).left(l).to_Maier()
    assert left2 == Maier_right, (N, l, left2, Maier_right)

    right2 = L.from_Maier(t).right(l).to_Maier()
    assert right2 == Maier_left, (N, l, right2, Maier_left)
