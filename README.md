# Modular Curves Framework

This repository contains framework for working with genus-0 modular curves, including all $X_0(N)$ curves and a few others.

## Included programs

- [test_relations.py](./test_relations.py) contains verification of the properties of the maps/equations, both symbolically (as rational maps) and on random values. To be ran with `sage -python -m pytest test_relations.py` (requires `pytest` installed in SageMath).
- [AlgebraicAttack.ipynb](./AlgebraicAttack.ipynb) implements a naive version of the meet-in-the-middle (GCD) algebraic attack against $2^n$-isogeny path problem.

## Quick API example

```python
> from parameters import Level4, Level16
> from setting import Setting

> t = Level16()
> t.right().right().right()

<Level2:D=256*r^8 - 256*r^4>

> t = Level4(6, S=Setting(p=19))
> t.is_supersingular()

True

# sample a 2-isogeny path of 5 elements
# following this parameter
> t.sample_fw(5)

[<Level4:A=3*i + 10>,
<Level4:A=10*i>,
<Level4:A=0>,
<Level4:A=6>,
<Level4:A=3*i + 10>]
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

- `t.left(l)`, `t.right(l)`
- `t.dual(l)`
- `t1.merge(t2, l)`
- `t.turn_tail_one()`, `t.turn_head_one()`
- `t.turn_tail_two()`, `t.turn_head_two()`