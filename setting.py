from sage.all import (
    proof, is_prime, GF, PolynomialRing, ZZ, QQ,
    sqrt, EllipticCurve, discrete_log,
    next_prime,
    legendre_symbol, fundamental_discriminant, hilbert_class_polynomial,
    NumberField, NumberFieldTower
)


proof.all(False)


class Setting:
    def __init__(self, p=None, pbits=None, p_complex=True):
        if p is None:
            p = 2**pbits
            while True:
                p = next_prime(p)
                if not p_complex or p % 4 == 3:
                    break

        self.p = ZZ(p)
        self.Fp = GF(self.p)

        if self.p % 4 == 3:
            # main supported case
            i = self.Fp['i'].gen()
            self.Fp2 = GF(self.p**2, name='i', modulus=i**2+1)
            self.i = self.Fp2([0, 1])
            assert str(self.i) == "i"

            #raw_A0 = self.Fp2(6)
            j0 = self.Fp2(1728)
        else:
            i = self.Fp(-1).sqrt()
            self.i = max(i, -i)
            self.Fp2 = self.Fp.extension(2)
            E = self.get_supersingular_curve()
            assert E.is_supersingular()
            j0 = E.j_invariant()  #.montgomery_model().a2()

        self.ONE = self.Fp2(1)

        self.g = self.Fp2.primitive_element()
        if (self.Fp2.order() - 1) % 3 == 0:
            self.ζ3 = self.Fp2(1).nth_root(3, all=True)[1]
            assert self.ζ3 != 1
        else:
            self.ζ3 = None

        if (self.Fp2.order() - 1) % 5 == 0:
            self.ζ5 = self.Fp2(1).nth_root(5, all=True)[1]
            assert self.ζ5 != 1
        else:
            self.ζ5 = None

        self.SQRT2 = self.Fp2(2).sqrt()
        try:
            self.CONST_TURN_TAIL25 = self.roots(lambda w: w**4 + 2*w**2 + self.ONE/5)[0]
            self.SQRT5 = 2 / (self.CONST_TURN_TAIL25**2 + 1)
        except (IndexError, ValueError):
            self.CONST_TURN_TAIL25 = NotImplemented
            self.SQRT5 = self.Fp2(5).sqrt()
        assert self.SQRT2**2 == 2
        assert self.SQRT5**2 == 5
        self.CONST_DUAL_X1_5 = self.Fp2(11)/2 + 5*self.SQRT5/2

        from parameters import Level1
        self.j0 = Level1(j0, S=self)
        self.D0 = self.j0.merge(self.j0.sample_fw(d=2), d=2)
        self.A0 = self.D0.merge(self.D0.sample_fw())
        self.a0 = self.A0.merge(self.A0.sample_fw())
        self.r0 = self.a0.merge(self.a0.sample_fw())

    def get_supersingular_curve(self):
        # https://eprint.iacr.org/2023/106
        # https://github.com/friends-of-quaternions/deuring
        p = self.p
        if p % 4 == 3:
            q = 1
            E = EllipticCurve(GF(p), [1, 0])
            # E.set_order(p + 1)
            EE = E.change_ring(self.Fp2)
            return EE
        else:
            for q in map(ZZ, range(3,p,4)):
                if not q.is_prime():
                    continue
                if legendre_symbol(-q, p) == -1:
                    break
            else:
                assert False  # should never happen

            H = hilbert_class_polynomial(fundamental_discriminant(-q))
            jinv = H.change_ring(GF(p)).any_root()

            if jinv == 0:
                E = EllipticCurve(self.Fp2, [0,1])
            else:
                a = 27 * jinv / (4 * (1728-jinv))
                E = EllipticCurve(self.Fp2, [a,-a])
            # E.set_order(p + 1)
            EE = E.change_ring(self.Fp2)
            return EE

    def Fp2vars(self, names, order=None):
        if order is None:
            return PolynomialRing(self.Fp2, names=names).gens()
        return PolynomialRing(self.Fp2, names=names, order=order).gens()

    def roots(self, f):
        if type(f).__name__ == "function":
            x, = self.Fp2vars("x")
            poly = f(x)
        else:
            poly = f
        return poly.roots(multiplicities=False)

    # SECTION: Montgomery 2-/3-/4-isogenies
    # See https://github.com/microsoft/SIKE-challenges/blob/main/montgomery_arithmetic.sage

    def iso2_push(self, x_ker, x):
        return x*(x_ker*x-1)/(x-x_ker)

    def iso3_push(self, x_ker, x):
        # x * (s*x-1)**2 / (x-s)**2
        # r * (s*r**2-1)**2 / (r**2-s)**2
        # r * s * (r**2-1/s)**2 / (r**2-s)**2)
        return x*(x_ker*x-1)**2/(x-x_ker)**2

    def iso2_A(self, x_ker):
        # A = 2 - 4*x**2
        # A-2=-4x**2
        # x**2 = (2-A)/4
        return 2*(1-2*x_ker**2)

    def iso4_A(self, x_ker):
        return 2*(1-2*x_ker**4)

    def iso3_A(self, x_ker, A):
        return (6-6*x_ker**2+A*x_ker)*x_ker

    # SECTION: modular polynomials
    # See https://math.mit.edu/~drew/ClassicalModPolys.html

    def MODULAR2(self, X, Y):
        return X**3 + Y**3 - 162000*(X**2 + Y**2 ) + 1488*(X**2*Y + X*Y**2) - X**2*Y**2 + 8748000000*(X + Y) + 40773375*X*Y - 157464000000000

    def MODULAR3(self, X, Y):
        res = 0
        res += (X**1 * Y**0 + Y**1 * X**0) * 1855425871872000000000
        res += (X**1 * Y**1) * -770845966336000000
        res += (X**2 * Y**0 + Y**2 * X**0) * 452984832000000
        res += (X**2 * Y**1 + Y**2 * X**1) * 8900222976000
        res += (X**2 * Y**2) * 2587918086
        res += (X**3 * Y**0 + Y**3 * X**0) * 36864000
        res += (X**3 * Y**1 + Y**3 * X**1) * -1069956
        res += (X**3 * Y**2 + Y**3 * X**2) * 2232
        res += (X**3 * Y**3) * -1
        res += (X**4 * Y**0 + Y**4 * X**0) * 1
        return res

    def MODULAR4(self, X, Y):
        return X**6 - X**5*Y**4 + 2976*X**5*Y**3 - 2533680*X**5*Y**2 + 561444609*X**5*Y - 8507430000*X**5 - X**4*Y**5 + 7440*X**4*Y**4 + 80967606480*X**4*Y**3 + 1425220456750080*X**4*Y**2 + 1194227244109980000*X**4*Y + 24125474716854750000*X**4 + 2976*X**3*Y**5 + 80967606480*X**3*Y**4 + 2729942049541120*X**3*Y**3 - 914362550706103200000*X**3*Y**2 + 12519806366846423598750000*X**3*Y - 22805180351548032195000000000*X**3 - 2533680*X**2*Y**5 + 1425220456750080*X**2*Y**4 - 914362550706103200000*X**2*Y**3 + 26402314839969410496000000*X**2*Y**2 + 188656639464998455284287109375*X**2*Y + 158010236947953767724187500000000*X**2 + 561444609*X*Y**5 + 1194227244109980000*X*Y**4 + 12519806366846423598750000*X*Y**3 + 188656639464998455284287109375*X*Y**2 - 94266583063223403127324218750000*X*Y - 364936327796757658404375000000000000*X + Y**6 - 8507430000*Y**5 + 24125474716854750000*Y**4 - 22805180351548032195000000000*Y**3 + 158010236947953767724187500000000*Y**2 - 364936327796757658404375000000000000*Y + 280949374722195372109640625000000000000

    # SECTION: create EllipticCurves

    def A_to_E(self, A):
        return EllipticCurve(self.Fp2, [0, A, 0, 1, 0])

    def lam_to_E(self, lam):
        return EllipticCurve(self.Fp2, [0, -lam-1, 0, lam, 0])


    # SECTION: research helpers

    def print_magma_system(self, eqs):
        for eq in eqs:
            vs = eq.parent().gens()

        print(f"""
F := GF({self.p});
Fi<ii> := PolynomialRing(F);
Fp2 := ext< F | ii^2 + 1>;
i := Fp2 ! ii;
P<{",".join(str(x) for x in vs)}> := PolynomialRing(Fp2, {len(vs)}, "grevlex");
I := ideal<P | """, end="")
        for eq in eqs:
            print(eq, ",", end="")
        print("0>;")
        print("B := GroebnerBasis(I); B")

    def to_QQ(self, v, lim=10**9, dlim=10**4, field=None, emax=9):
        if field is None:
            field = self.Fp2
        v = field(v)

        if emax <= 1:
            if int(v) < lim:
                return v
            if int(-v) < lim:
                return "-%d" % int(-v)
            if int(1/v) < lim:
                return "Fp2(1)/%d" % int(1/v)
            if int(-1/v) < lim:
                return "-Fp2(1)/%d" % int(-1/v)
            for d in range(1, dlim):
                if int(v*d) < lim:
                    return "Fp2(%d)/%d" % (int(v*d), d)
                if int(-v*d) < lim:
                    return "-Fp2(%d)/%d" % (int(-v*d), d)
            for d in range(1, dlim):
                if int(d/v) < lim:
                    return "Fp2(%d)/%d" % (d, int(d/v))
                if int(-d/v) < lim:
                    return "-Fp2(%d)/%d" % (-d, int(d/v))
            raise
        else:
            for e in range(1, emax+1):
                try:
                    ret = self.to_QQ(v**e, lim=lim, dlim=dlim, field=field, emax=0)
                except:
                    continue
                if e == 1:
                    return ret
                return f"root{e}({ret})"
            raise

    def poly_to_QQ(self, poly, **opts):
        opts.setdefault("lim", int(sqrt(self.p)))
        opts.setdefault("dlim", 258)
        res = []
        target = poly
        for c, mono in target:
            cc = self.to_QQ(c, **opts)
            res.append("%s*%s" % (cc, mono))
        return (" " + " + ".join(res)).replace("^", "**").replace("+ -", "- ").replace(" 1*", " ").strip()

    def poly_to_fraction(self, f, var):
        # var = top / bot
        # var * bot - top = 0
        top = -f.subs({var: 0})
        bot = (f+top).subs({var: 1})
        assert f == var*bot - top
        return top, bot


class RationalSetting(Setting):
    def __init__(self):
        self.Fp2 = self.Fp = QQ
        self.ONE = self.Fp2(1)


class CubicSetting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberFieldTower([x**2+x+1], names="z3")
        self.ζ3 = self.Fp2.gen()
        self.ONE = self.Fp2(1)


class ComplexSetting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberField(x**2+1, name="i")
        self.i = self.Fp2.gen()
        self.ONE = self.Fp2(1)

class Sqrt5Setting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberField(x**2-5, name="s5")
        self.SQRT5 = self.Fp2.gen()
        self.ONE = self.Fp2(1)

class Sqrt5Zeta5Setting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberFieldTower([
            x**4+x**3+x**2+x+1, # zeta5
            x**2+x+1, # zeta3
            # x**2-5, # sqrt(5)
        ], names="z5,z3") #,s5")
        #self.i, self.ζ5, self.ζ3, self.SQRT5 = self.Fp2.gens()
        self.ζ5, self.ζ3 = self.Fp2.gens()

        self.SQRT5 = 2*self.ζ5**3 + 2*self.ζ5**2 + 1
        self.ONE = self.Fp2(1)

        assert self.ζ3**3 == 1
        assert self.ζ5**5 == 1
        assert self.SQRT5**2 == 5


class ComplexCubicSetting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberFieldTower([x**2+1, x**2+x+1], names="i,z3")
        self.i, self.ζ3 = self.Fp2.gens()
        self.ONE = self.Fp2(1)

class ComplexCubicSqrt5Setting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberFieldTower([x**2+1, x**2+x+1,x**2-5], names="i,z3,s5")
        self.i, self.ζ3, self.SQRT5 = self.Fp2.gens()
        self.ONE = self.Fp2(1)

class ComplexCubicSqrt5Zeta5Setting(RationalSetting):
    def __init__(self):
        x = QQ['x'].gen()
        self.Fp2 = NumberFieldTower([
            x**4+x**3+x**2+x+1, # zeta5
            x**2+1, # sqrt(-1)
            x**2+x+1, # zeta3
            # x**2-5, # sqrt(5)
        ], names="z5,i,z3") #,s5")
        #self.i, self.ζ5, self.ζ3, self.SQRT5 = self.Fp2.gens()
        self.ζ5, self.i, self.ζ3 = self.Fp2.gens()

        self.SQRT5 = 2*self.ζ5**3 + 2*self.ζ5**2 + 1
        self.ONE = self.Fp2(1)

        assert self.i**2 == -1
        assert self.ζ3**3 == 1
        assert self.ζ5**5 == 1
        assert self.SQRT5**2 == 5
