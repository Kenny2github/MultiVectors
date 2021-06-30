Welcome to MultiVectors documentation!

- [Concepts](#concepts)
  - [Bases and geometric products](#bases-and-geometric-products)
  - [Blades and multivectors](#blades-and-multivectors)
  - [The choose operator, inner (dot) and outer (wedge) products](#the-choose-operator-inner-dot-and-outer-wedge-products)
  - [Euler's formula applied to multivectors](#eulers-formula-applied-to-multivectors)
- [Applied](#applied)
  - [Blade](#blade)
  - [MultiVector](#multivector)

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
* For reasons that are beyond my power to explain, the rotation of a multivector *V* by *θ* through the plane *B* is e\*\*(-*θB*/2) * V * e\*\*(*θB*/2)

## Applied
All of the above concepts are applied in this library.

### Blade
* `multivectors.Blade(*bases, scalar=a)` represents a blade with **0-indexed bases** `bases` multiplied by a real scalar `a`. For example, `Blade(0, 1, 3)` represents the basis volume *x̂ŷŵ*, with 0 meaning *x̂*.
* Basis names can be *swizzled* on the `Blade` class itself: the above could have been done with `Blade.xyw`. For basis vectors beyond *ŵ*, use *ê*ₙ: `Blade.e1e3e4e5` is a basis 4-vector in 5D space. `Blade._` is a 0-blade: a scalar, but with `Blade` type.
* Basis indices can also be used: `Blade[:4]` = `Blade[0, 1, 2, 3]` = `Blade(0, 1, 2, 3)` = `Blade.xyzw`
* To take it one step further, basis names can be swizzled on the module itself: `multivector.xyz` returns `multivectors.Blade.xyz`. This also works when importing: `from multivectors import x, y, z, xy, xz, yz` will work just fine.
* The rules of arithmetic with blades as described above apply:
```python
>>> from multivectors import x, y, z
>>> x * 2 + x * 3
5.0 * Blade.x
>>> x * 5 * y
5.0 * Blade.xy
>>> (x * y) * z == x * (y * z)
True

```
* You can query the grade of a blade: `Blade.xy.grade` is 2.

### MultiVector
* `multivectors.MultiVector.from_terms(*terms)` represents a **sum of `terms`**. You normally should not be constructing this class; it is created from summation involving `Blade`s.
* Basis names can be swizzled on class **instances** to get the coefficient of that basis:
```python
>>> from multivectors import x, y
>>> (x + 2*y).x
1.0

```
* Basis indices can also be used:
```python
>>> from multivectors import xy, yz
>>> (xy + 2*yz)[0, 1]
1.0

```
* Choosing by grade is supported:
```python
>>> from multivectors import x, y, xy, yz
>>> V = 1 + 2*x + 3*y + 4*xy + 5*yz
>>> V % 0
1.0
>>> V % 1
(2.0 * Blade.x + 3.0 * Blade.y)
>>> V % 2
(4.0 * Blade.xy + 5.0 * Blade.yz)

```
* The rules of arithmetic with multivectors as described above apply:
```python
>>> from multivectors import x, y, z
>>> x * (y + z)
(1.0 * Blade.xy + 1.0 * Blade.xz)
>>> (y + z) * (2 * x)
(-2.0 * Blade.xy + -2.0 * Blade.xz)
>>> (1*x + 2*y) * (3*x + 4*y)
(11.0 + -2.0 * Blade.xy)

```
* The extra products apply too:
```python
>>> from multivectors import x, y, z
>>> 1*3 + 2*4
11
>>> (1*x + 2*y) @ (3*x + 4*y)
11.0
>>> (1*5 - 2*4, 1*6 - 3*4, 2*6 - 3*5)
(-3, -6, -3)
>>> (1*x + 2*y + 3*z) ^ (4*x + 5*y + 6*z)
(-3.0 * Blade.xy + -6.0 * Blade.xz + -3.0 * Blade.yz)

```
* A convenience method is provided to rotate multivectors:
```python
>>> from math import radians
>>> from multivectors import x, y, z, xz
>>> round((3*x + 2*y + 4*z).rotate(radians(90), xz), 2)
(-4.0 * Blade.x + 2.0 * Blade.y + 3.0 * Blade.z)

```