import pytest
import warnings

from sage.all import is_prime, divisors, gcd

from setting import *
from parameters import *

def WARN(*args):
    warnings.warn(" ".join(map(str, args)))

LEVELS = [
    Level1,

    Level2,
    Level4,
    Level8,
    Level16,

    Level3,
    Level6,
    Level9,
    Level10,
    Level12,
    Level18,

    Gamma2,
    Gamma2_Level2,
    Gamma3,
    Gamma3_Level3,

    Level5,
    Level5_X1,
    Level25,
    Level25_X5,

    Level7,
    Level13,
]
assert set(LEVELS) == set(Parameter.__subclasses__())
LEVEL_PRIME_ELL_PAIRS = []
LEVEL_ANYDIV_ELL_PAIRS = []
LEVEL_PRIME_ELL_ELL2_PAIRS = []
LEVEL_ANYDIV_ELL_ELL2_PAIRS = []
for level in LEVELS:
    for ell in divisors(level.LEVEL):
        if ell > 1:
            LEVEL_ANYDIV_ELL_PAIRS.append((level, ell))
            if is_prime(ell):
                LEVEL_PRIME_ELL_PAIRS.append((level, ell))
            for ell2 in divisors(level.LEVEL):
                if ell < ell2 and gcd(ell, ell2) == 1:
                    LEVEL_ANYDIV_ELL_ELL2_PAIRS.append((level, ell, ell2))
                    if is_prime(ell) and is_prime(ell2):
                        LEVEL_PRIME_ELL_ELL2_PAIRS.append((level, ell, ell2))


@pytest.mark.parametrize('level,l', LEVEL_ANYDIV_ELL_PAIRS)
def test_dual_is_involution(level, l):
    t = level()
    # dual is involution
    assert t.dual(l).dual(l) == t
    assert t.dual().dual() == t


@pytest.mark.parametrize('level', LEVELS)
def test_dual_reverses_J(level):
    t = level()
    # dual acts on J(t) by reversal
    assert t.j_seq() == t.dual().j_seq()[::-1]


@pytest.mark.parametrize('level,l', LEVEL_ANYDIV_ELL_PAIRS)
def test_left_right_anticommute_with_dual(level, l):
    t = level()
    # L,R anticommute with respective duals
    assert t.dual(l).left(l) == t.right(l).dual(l, ignore_extra_level=True)
    assert t.dual(l).left(l).dual(l, ignore_extra_level=True) == t.right(l)
    assert t.dual(l).right(l).dual(l, ignore_extra_level=True) == t.left(l)


@pytest.mark.parametrize('level,l,ll', LEVEL_ANYDIV_ELL_ELL2_PAIRS)
def test_left_right_anticommute_coprime(level, l, ll):
    t = level()
    assert t.left(l).right(ll) == t.right(ll).left(l)


@pytest.mark.parametrize('level,l', LEVEL_ANYDIV_ELL_PAIRS)
def test_phi_left_right_is_0(level, l):
    """
    Phi_{N/l,l}( L_{N,l}(t), R_{N,l}(t) ) = 0
    """
    t = level()
    if hasattr(t.left(l), "raw_phi_eq%d" % l):
        assert t.left(l).phi(t.right(l), l) == 0
    else:
        # for composite l we do not derive Phi usually
        if is_prime(l):
            WARN(type(t.left(l)).__name__, "no phi eq for", l)


@pytest.mark.parametrize('level,l', LEVEL_PRIME_ELL_PAIRS)
def test_left_right_merge(level, l):
    """
    M(L,R)=Id
    """
    t = level()
    if level == Level5_X1:
        WARN("skipping X_1(5) LR-Merge test - no merge map from X(1)")
        return
    assert t == t.left(l).merge(t.right(l), l=l)
