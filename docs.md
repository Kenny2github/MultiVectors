Welcome to MultiVectors documentation!

- [Concepts](#concepts)
  - [Bases and geometric products](#bases-and-geometric-products)
  - [Blades and multivectors](#blades-and-multivectors)
  - [The choose operator, inner (dot) and outer (wedge) products](#the-choose-operator-inner-dot-and-outer-wedge-products)
  - [Euler's formula applied to multivectors](#eulers-formula-applied-to-multivectors)
- [Applied](#applied)
  - [Blades](#blades)
  - [MultiVector](#multivector)
  - [Scalar factory](#scalar-factory)

## Concepts
Here are some concepts to bear in mind. This is a *very* brief introduction to geometric algebra; a longer one is [here](https://www.youtube.com/watch?v=60z_hpEAtD8).

### Bases and geometric products
* Every dimension of space comes with a *basis vector*: an arrow of length 1 unit pointed towards the positive end of the axis.
  * In our 3-dimensional world, there are the basis vectors *x̂*, *ŷ*, and *ẑ*, which are pointed towards the positive ends of the *x*, *y*, and *z* axes respectively.
  * The fourth dimensional basis vector is *ŵ*. In higher dimensions, usually all bases are numbered instead of lettered: the fifth dimensional basis vectors are *ê₁*, *ê₂*, *ê₃*, *ê₄*, and *ê₅*.
* The *geometric product* of two basis vectors is their simple multiplication - not the dot or cross product! The geometric product of *x̂* and *ŷ* is simply *x̂ŷ*.
  * The geometric product of two basis vectors is a *basis plane*. *x̂ŷ* is the basis plane of the *x-y* plane. The other basis planes are *ŷẑ* and *x̂ẑ*.
  * The geometric product of three basis vectors is a *basis volume*. *x̂ŷẑ* is the basis volume of 3D space, which only has one basis volume, but also one of the four basis volumes of 4D space.
* The geometric product of a basis vector with itself is 1. That is, *x̂x̂* = *x̂*² = *ŷŷ* = *ŷ*² = *ẑẑ* = *ẑ*² = 1
* The geometric product of different basis vectors *anticommutes*: *x̂ŷ* = -*ŷx̂* and *x̂ŷẑ* = -*x̂ẑŷ* = *ẑx̂ŷ* = -*ẑŷx̂*

### Blades and multivectors
* A *blade* is a *scaled basis*: a *scalar* (regular real number) multiplied by a *basis*. For example, 3*x̂ŷ* is a blade. Note that this means all bases are blades scaled by 1.
  * A *k-blade* is a blade of *grade k*: the geometric product of a scalar and *k* different basis vectors. 3*x̂ŷ* has grade 2; it is a 2-blade.
  * Scalars are 0-blades - blades consisting of *no* basis vectors.
* A *multivector* is a sum of multiple blades. For example, 1 + 2*x̂* - 3*ŷẑ* is a multivector.
  * The sum of multiple (and only) 1-blades is usually called a simple *vector*. For example, 3*x̂* + 2*ŷ* is a vector.
  * The sum of multiple (and only) 2-blades is a *bivector*. Basis planes are also known as *basis bivectors*. For example, 3*x̂ŷ* is a bivector.
* The rules of linearity, associativity and distributivity in multiplication apply, as long as order of arguments is maintained:
  * (*x̂*)(*aŷ*) = *ax̂ŷ* (linearity, for scalar *a*)
  * (*x̂ŷ*)*ẑ* = *x̂*(*ŷẑ*) (associativity)
  * *x̂*(*ŷ* + *ẑ*) = *x̂ŷ* + *x̂ẑ* (distributivity)
  * (*ŷ* + *ẑ*)(*ax̂*) = *a*(*ŷ* + *ẑ*)(*x̂*) (linearity) = a(*ŷx̂* + *ẑx̂*) (distributivity) = a(-*x̂ŷ* - *x̂ẑ*) (anticommutativity)
* However, some things which require commutativity break down, such as the binomial theorem.

### The choose operator, inner (dot) and outer (wedge) products
* ⟨*V*⟩ₙ *chooses* all *n*-blades from the multivector *V*. For example, if *V* = 1 + 2*x̂* + 3*ŷ* + 4*x̂ŷ* + 5*ŷẑ*, then ⟨*V*⟩₀ = 1 and ⟨*V*⟩₁ = 2*x̂* + 3*ŷ* and ⟨*V*⟩₂ = 4*x̂ŷ* + 5*ŷẑ*
* *U* · *V* = ⟨*UV*⟩ₙ where *U* is of grade *r*, *V* is of grade *s*, and *n* = |*r - s*|. This is the *inner* or *dot product*.
  * The dot product associates and distributes the same way the geometric product does.
  * From this, for arbitrary vectors *ax̂* + *bŷ* and *cx̂* + *dŷ*, we recover the typical meaning of the dot product:<br/>![Derivation of normal vector dot product](https://cdn.discordapp.com/attachments/417244106876780544/859663520747356170/dot_product.png)
* *U* ∧ *V* = ⟨*UV*⟩ₙ where *U* is of grade *r*, *V* is of grade *s*, and *n* = *r + s*. This is the *outer* or *wedge product*.
  * The outer product associates and distributes the same way the geometric product does.
  * From this, for arbitrary vectors *ax̂* + *bŷ* + *cẑ* and *dx̂* + *eŷ* + *fẑ*, we recover something that looks very much like a cross product:<br/>![Derivation of normal vector cross product](https://cdn.discordapp.com/attachments/417244106876780544/859663658762108968/wedge_product.png)

### Euler's formula applied to multivectors
* *e* to the power (*θB*) = cos(*θ*) + *B* sin(*θ*) where θ is a scalar in radians and *B* is a basis multivector.
  * **Note**: This formula only works when *B*² = -1 (in the same way as the imaginary unit *i*).
  * For purposes of interest, the more general formula for *e* raised to a multivector power is the Taylor series:<br/>![Taylor series of exp(V)](https://cdn.discordapp.com/attachments/417244106876780544/905703412488896522/unknown.png)
* For reasons that are beyond my power to explain, the rotation of a multivector *V* by *θ* through the plane *B* is e\*\*(-*θB*/2) * V * e\*\*(*θB*/2)

## Applied
All of the above concepts are applied in this library.

### Blades
* `multivectors.x`, `y`, `z`, `w` represent the first four basis vectors. Bases past 4D use *ê*ₙ syntax: `multivectors.e5` is the 5th basis vector.
* Basis names can be *swizzled* on the module:
```python
>>> from multivectors import x, y, z, xyz
>>> x * y * z == xyz
True

```
* `multivectors._` is a 0-blade: a scalar, but with `MultiVector` type. It is differentiated by normal scalars when `repr()`ed by surrounding parentheses:
```python
>>> from multivectors import _
>>> 2 * _
(2.0)

```
* The rules of arithmetic with blades as described above apply:
```python
>>> from multivectors import x, y, z
>>> x * 2 + x * 3
(5.0 * x)
>>> x * 5 * y
(5.0 * x*y)
>>> (x * y) * z == x * (y * z)
True

```
* You can query the grade of a blade (only if it is an actual blade!):
```python
>>> from multivectors import xy, yzw
>>> xy.grade
2
>>> (2 * yzw).grade
3

```

### MultiVector
* A multivector is a sum of zero or more blades. All scalars and blades are multivectors, but not all multivectors are scalars or blades.
* Basis names can be swizzled on class **instances** to get the coefficient of that basis:
```python
>>> from multivectors import x, y
>>> (x + 2*y).x
1.0

```
* Basis indices can also be used:
```python
>>> from multivectors import xy, yz
>>> (xy + 2*yz) % (0, 1)
1.0

```
* Choosing by grade is supported:
```python
>>> from multivectors import x, y, xy, yz
>>> V = 1 + 2*x + 3*y + 4*xy + 5*yz
>>> V[0]
(1.0)
>>> V[1]
(2.0 * x + 3.0 * y)
>>> V[2]
(4.0 * x*y + 5.0 * y*z)

```
* Getting the grade of a multivector will only work if if the multivector is a blade. Otherwise, you get `None`:
```python
>>> from multivectors import w, xy, yz
>>> w.grade
1
>>> xy.grade
2
>>> yz.grade
2
>>> (w + xy).grade
>>> # returns None
>>> (xy + yz).grade
>>> # even multivectors consisting of blades with the same grade don't work
>>> # the multivector must be an actual blade, with only one basis sequence

```
* The rules of arithmetic with multivectors as described above apply:
```python
>>> from multivectors import x, y, z
>>> x * (y + z)
(1.0 * x*y + 1.0 * x*z)
>>> (y + z) * (2 * x)
(-2.0 * x*y + -2.0 * x*z)
>>> (1*x + 2*y) * (3*x + 4*y)
(11.0 + -2.0 * x*y)

```
* The extra products apply too:
```python
>>> from multivectors import x, y, z
>>> 1*3 + 2*4
11
>>> (1*x + 2*y) @ (3*x + 4*y)
(11.0)
>>> (1*5 - 2*4, 1*6 - 3*4, 2*6 - 3*5)
(-3, -6, -3)
>>> (1*x + 2*y + 3*z) ^ (4*x + 5*y + 6*z)
(-3.0 * x*y + -6.0 * x*z + -3.0 * y*z)

```
* A convenience method is provided to rotate multivectors:
```python
>>> from math import radians
>>> from multivectors import x, y, z, xz
>>> round((3*x + 2*y + 4*z).rotate(radians(90), xz), 2)
(-4.0 * x + 2.0 * y + 3.0 * z)

```

### Scalar factory
You can set a scalar factory to control the type of number used for scalars. This can be, for example:
* [`decimal.Decimal`](https://docs.python.org/3/library/decimal.html) to get arbitrary decimal precision in multivectors
* [`fractions.Fraction`](https://docs.python.org/3/library/fractions.html) to limit calculations to rationals
* A class of your own, to control other aspects of scalar arithmetic

You will need at least two things:
* `instancecheck`: The class object(s) to limit arithmetic to.
  * This can be a tuple of classes, and should include [`numbers.Real`](https://docs.python.org/3/library/numbers.html#numbers.Real) if you want to retain regular literal arithmetic. If you only include the custom class, things like `2 * x` will fail because `int` is not a subclass of your class.
* `factory`: The factory function that creates instances of your class. This can just be the class object again, but note that the function/constructor must implement the following:
  * `factory()` returns the 0 value of the class. This is used in defaults.
  * `factory(f)` where `f` is a `float` value. This is used to cast floats to the class.
  * `factory(s)` where `s` is a `str` value. This is used for non-zero constants like `factory('1')`.
* Optionally, `annotation`: The annotation used in type hinting.

You must call `multivectors.set_scalar_factory()` **before importing anything else, including swizzled names**:
```python
>>> from decimal import Decimal
>>> from numbers import Real
>>> from multivectors import set_scalar_factory
>>> set_scalar_factory((Decimal, Real), Decimal)
>>> from multivectors import x, y
>>> x + y
(Decimal('1') * x + Decimal('1') * y)

```
<div style="display: none">
Resetting scalar factory for doctest:
```python
>>> set_scalar_factory(Real, float, float)

```
</div>