# Modular Curves Framework

This repository contains framework for working with genus-0 modular curves, including all $X_0(N)$ curves and a few others.

## Included programs

- [test_relations.py](./test_relations.py) contains verification of the properties of the maps/equations, both symbolically (as rational maps) and on random values. To be ran with `sage -python -m pytest test_relations.py` (requires `pytest` installed in SageMath).
- [AlgebraicAttack.ipynb](./AlgebraicAttack.ipynb) implements a naive version of the meet-in-the-middle (GCD) algebraic attack against $2^n$-isogeny path problem.
- [Examples.ipynb](./Examples.ipynb) contains some usage examples in the form of notebook.

## Quick API example

```python
> from parameters import Level1, Level4, Level8, Level16
> from setting import Setting
> t = Level16()
> t.left(), t.right()

(<Level8:a=(1/2*r^2 + 1/2)/r>, <Level8:a=r^2>)

> t.left().phi(t.right())

0

> x, y = QQ['x,y'].gens()
> Level16(x).phi(Level16(y))

-2*x^2*y + y^2 + 1

> t.left().merge(t.right())

<Level16:r=r>

> Level8(x).merge(Level8(y), check=False)
<Level16:r=(y + 1)/(2*x)>

> t = Level4(6, S=Setting(p=19))
> t.is_supersingular()

True

> # sample a 2-isogeny path of 5 elements
> # following this parameter
> t.sample_fw(5, skip=0)

[<Level4:A=6>,
<Level4:A=16*i + 10>,
<Level4:A=9*i>,
<Level4:A=13>,
<Level4:A=16*i + 10>]

> # propagate 2- and 3- isogeny from j-invariants
> j1 = Setting(p=1051).j0.sample_fw(1, l=2, skip=10) # do some walk to get a random one
> j2 = j1.sample_fw(l=2)
> j3 = j1.sample_fw(l=3)
> assert j1.phi(j2, l=2) == 0
> assert j1.phi(j3, l=3) == 0
> j1, j2, j3

(<Level1:j=403*i + 524>, <Level1:j=798*i + 457>, <Level1:j=443*i + 842>)

> # note: reverse order
> # because MergeBar requires
> # equality of rights
> D12 = j2.merge(j1, l=2)
> D13 = j3.merge(j1, l=3)
> D6 = D12.merge(D13)
> D6

<Level6:t=480*i + 481>

> j6 = D6.left().left()
> assert j6.phi(j2, l=3) == 0
> assert j6.phi(j3, l=2) == 0
> j6

<Level1:j=390*i + 6>
```

## Code structure

[setting.py](./setting.py) contains `Setting` class which corresponds to the choice of the field $GF(p^2)$ and implements many general useful functions. There are also other `Settings` classes which work over certain number fields; their main purpose is to provide algebraic numbers (e.g. roots of unity) for symbolic computations over extensions of $\mathbb{Q}$. For example, `ComplexSetting` provides $\sqrt{-1}$, `ComplexCubicSqrt5Zeta5Setting` provides roots of unity $\zeta_3,\zeta_4,\zeta_5$ and $\sqrt{5}=2\zeta_5^3 + 2\zeta_5^2 + 1$.

[parameters.py](./parameters.py) is the main file containing classes for parameters of modular curves, together with implementations of the corresponding maps and equations. Parameter $t_N \in X_0(N)$ is described by an instance of the `Level{N}` class (e.g., `Level2`, `Level18`, etc.).

### List of modular curves

The full list of currently supported parameters is:

Binary hierarchy:
- `Level1` ( $X_0(1)$ )
- `Level2` ( $X_0(2)$ )
- `Level4` ( $X_0(4)$ )
- `Level8` ( $X_0(8)$ )
- `Level16` ( $X_0(16)$ )

Multiples of 3:
- `Level3` ( $X_0(3)$ )
- `Level6` ( $X_0(6)$ )
- `Level9` ( $X_0(9)$ )
- `Level12` ( $X_0(12)$ )
- `Level18` ( $X_0(18)$ )

Non-split Cartan subgroups (inside towers $X_0(12)$ and $X_0(18)$ ):
- `Gamma3` ( $\gamma_3 \in X_{\mathrm{ns}}(2)$, $\gamma_3=\sqrt{j-1728}$ )
- `Gamma3_Level3` ( $Z=X_0(3)\otimes X_{\mathrm{ns}}(2)$ )
- `Gamma2` ( $\gamma_2 \in X_{\mathrm{ns}}^+(3)$, $\gamma_2=\sqrt[3]{j}$ )
- `Gamma2_Level2` ( $Y=X_0(2)\otimes X_{\mathrm{ns}}^+(3)$ )

Multiples of 5:
- `Level5` ( $X_0(5)$ )
- `Level10` ( $X_0(10)$ )
- `Level25` ( $X_0(25)$ )
- `Level5_X1` ( $X_1(5)$ )
- `Level25_X5` ( $X(5)$ )

Remaining $X_0(N)$ of genus 0:
- `Level7` ( $X_0(7)$ )
- `Level13` ( $X_0(13)$ )


### Instantiating parameters

The parameter can be instantiated in the following ways:

- Create symbolic parameter over $\mathbb{Q}$:\
  `t = Level16()`
- Create symbolic parameter over a given number field:\
  `t = Level16(S=ComplexSetting())`
- Create parameter over a given finite field $\mathbb{F}_{p^2}$:\
  `t = Level0(1728, S=Setting(p=19))`

The value of the parameter (symbolic or a field element) can be retried as `t.value`.

### List of main methods

Let $t,t_1,t_2 \in X_0(N)$ (instances of `Level{N}`).

- `t.left(l)`, `t.right(l)`: covers $\mathbf{L}\_{N,\ell}(t),\mathbf{R}\_{N,\ell}(t)$;
- `t.dual()`: full duality (Fricke) $\mathbf{w}_{N,N}(t)$;
- `t.dual(l)`: $\ell$-duality (Atkin-Lehner) $\mathbf{w}_{N,\ell}(t)$;
- `t1.merge(t2, l)`: inverse of covers $\mathbf{M}_{N,\ell}(t)$;
- `t.turn_tail_one(l)`, `t.turn_head_one(l)`: single-edge head/tail rotations $\mathbf{T}\_{N,\ell}(t),\mathbf{H}\_{N,\ell}(t)$
- `t.turn_tail_two(l)`, `t.turn_head_two(l)`: double-edge head/tail rotations $\mathbf{T}\_{N,\ell^2}(t),\mathbf{H}\_{N,\ell^2}(t)$