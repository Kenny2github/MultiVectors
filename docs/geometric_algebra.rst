.. _geo-intro:

Quick Intro to Geometric Algebra
================================

Here are some concepts to bear in mind. This is a *very* brief introduction to geometric algebra; a longer one is `here <https://www.youtube.com/watch?v=60z_hpEAtD8>`_.

.. default-role:: math

Bases and geometric products
----------------------------

* Every dimension of space comes with a *basis vector*: an arrow of length 1 unit pointed towards the positive end of the axis.

  * In our 3-dimensional world, there are the basis vectors `\hat x`, `\hat y`, and `\hat z`, which are pointed towards the positive ends of the `x`, `y`, and `z` axes respectively.
  * The fourth dimensional basis vector is `\hat w`. In higher dimensions, usually all bases are numbered instead of lettered: the fifth dimensional basis vectors are `\hat e_1`, `\hat e_2`, `\hat e_3`, `\hat e_4`, and `\hat e_5`.

* The *geometric product* of two basis vectors is their simple multiplication - not the dot or cross product! The geometric product of `\hat x` and `\hat y` is simply `\hat x \hat y`.

  * The geometric product of two basis vectors is a *basis plane*. `\hat x \hat y` is the basis plane of the `xy` plane. The other basis planes are `\hat y \hat z` and `\hat x \hat z`.
  * The geometric product of three basis vectors is a *basis volume*. `\hat x \hat y \hat z` is the basis volume of 3D space, which only has one basis volume, but also one of the four basis volumes of 4D space.

* The geometric product of a basis vector with itself is 1. That is, `\hat x \hat x = \hat x^2 = \hat y \hat y = \hat y^2 = \hat z \hat z = \hat z^2 = 1`.
* The geometric product of different basis vectors *anticommutes*: `\hat x \hat y = - \hat y \hat x` and `\hat x \hat y \hat z = - \hat x \hat z \hat y = \hat z \hat x \hat y = - \hat z \hat y \hat x`.

Blades and multivectors
-----------------------

* A *blade* is a *scaled basis*: a *scalar* (regular real number) multiplied by a *basis*. For example, `3 \hat x \hat y` is a blade. Note that this means all bases are blades scaled by 1.

  * A `k`-*blade* is a blade of *grade* `k`: the geometric product of a scalar and `k` different basis vectors. `3 \hat x \hat y` has grade `2`; it is a `2`-blade.
  * Scalars are `0`-blades - blades consisting of *no* basis vectors.

* A *multivector* is a sum of multiple blades. For example, `1 + 2 \hat x - 3 \hat y \hat z` is a multivector.

  * The sum of multiple (and only) 1-blades is usually called a simple *vector*. For example, `3 \hat x + 2 \hat y` is a vector.
  * The sum of multiple (and only) 2-blades is a *bivector*. Basis planes are also known as *basis bivectors*. For example, `3 \hat x \hat y` is a bivector.

* The rules of linearity, associativity and distributivity in multiplication apply, as long as order of arguments is maintained:

  * `(\hat x)(a \hat y) = a \hat x \hat y` (linearity, for scalar `a`)
  * `(\hat x \hat y)(\hat z) = \hat x (\hat y \hat z)` (associativity)
  * `\hat x (\hat y + \hat z) = \hat x \hat y + \hat x \hat z` (distributivity)
  * `(\hat y + \hat z)(a \hat x) = a (\hat y + \hat z)(\hat x)` (linearity) `= a (\hat y \hat x + \hat z \hat x)` (distributivity) `= a (- \hat x \hat y - \hat x \hat z)` (anticommutativity) `= -a (\hat x \hat y + \hat x \hat z)` (converse of distributivity)

* However, some things which require commutativity break down, such as the binomial theorem.

The choose operator, inner (dot) and outer (wedge) products
-----------------------------------------------------------

* `\langle V \rangle_n` *chooses* all `n`-blades from the multivector `V`. For example, if `V = 1 + 2 \hat x + 3 \hat y + 4 \hat x \hat y + 5 \hat y \hat z`, then `\langle V \rangle_0 = 1` and `\langle V \rangle_1 = 2 \hat x + 3 \hat y` and `\langle V \rangle_2 = 4 \hat x \hat y + 5 \hat y \hat z`.
* `U \cdot V = \langle UV \rangle_n` where `U` is of grade `r`, `V` is of grade `s`, and `n = |r - s|`. This is the *inner* or *dot product*.

  * The dot product associates and distributes the same way the geometric product does.
  * From this, for arbitrary vectors `a \hat x + b \hat y` and `c \hat x + d \hat y`, we recover the typical meaning of the dot product:

    .. math::

        & (a \hat x + b \hat y) \cdot (c \hat x + d \hat y) \\
        &= ac (\hat x \cdot \hat x) + ad (\hat x \cdot \hat y) + bc (\hat y \cdot \hat x) + bd (\hat y \cdot \hat y) \\
        &= ac \langle \hat x \hat x \rangle_0 + ad \langle \hat x \hat y \rangle_0 + bc \langle \hat y \hat x \rangle_0 + bd \langle \hat y \hat y \rangle_0 \\
        &= ac \langle 1 \rangle_0 + ad (0) + bc (0) + bd \langle 1 \rangle_0 \\
        &\text{(because }\hat x \hat y\text{ and }\hat y \hat z\text{ have no part with grade }0\text{)} \\
        &= ac + bd

* `U \wedge V = \langle UV \rangle_n` where `U` is of grade `r`, `V` is of grade `s`, and `n = r + s`. This is the *outer* or *wedge product*.

  * The outer product associates and distributes the same way the geometric product does.
  * From this, for arbitrary vectors `a \hat x + b \hat y + c \hat z` and `d \hat x + e \hat y + f \hat z`, we recover something that looks very much like a cross product:

    .. math::

        &(a \hat x + b \hat y + c \hat z) \wedge (d \hat x + e \hat y + f \hat z) \\
        &= ad (\hat x \wedge \hat x) + bd (\hat y \wedge \hat x) + cd (\hat z \wedge \hat x) \\
        &\quad + ae (\hat x \wedge \hat y) + be (\hat y \wedge \hat y) + ce (\hat z \wedge \hat y) \\
        &\quad + af (\hat x \wedge \hat z) + bf (\hat y \wedge \hat z) + cf (\hat z \wedge \hat z) \\
        &= ad \langle \hat x \hat x \rangle_2 + be \langle \hat y \hat y \rangle_2 + cf \langle \hat z \hat z \rangle_2 \\
        &\quad + (ae - bd) \langle \hat x \hat y \rangle_2 + (af - cd) \langle \hat x \hat z \rangle_2 + (bf - ce) \langle \hat y \hat z \rangle_2 \\
        &= 0 + 0 + 0 + \begin{vmatrix} a & b \\ d & e \end{vmatrix} (\hat x \hat y) + \begin{vmatrix} a & c \\ d & f \end{vmatrix} (\hat x \hat z) + \begin{vmatrix} b & c \\ e & f \end{vmatrix} (\hat y \hat z) \\
        &\text{(because }\hat x \hat x\text{ etc.} = 1\text{, which has no grade }2\text{ part)} \\
        &= \begin{vmatrix} \hat y \hat z & \hat x \hat z & \hat x \hat y \\ a & b & c \\ d & e & f \end{vmatrix}

Euler's formula applied to multivectors
---------------------------------------

* `e^{\theta B} = \cos \theta + B \sin \theta` where `\theta` is a scalar in radians and `B` is a basis multivector.

  .. note:: This formula only works when `B^2 = -1` (in the same way as the imaginary unit `i`), such as `(\hat x \hat y)^2 = \hat x \hat y \hat x \hat y = - \hat x \hat y \hat y \hat x = - \hat x \hat x = -1`.

  * For purposes of interest, the more general formula for `e` raised to a multivector power is the Taylor series:

    .. math::

        e^V = \exp(V) = \sum_{n=0}^{\infty} \frac{V^n}{n!}

* For reasons that are beyond my power to explain, the rotation of a multivector `V` by `\theta` through the plane `B` is `e^{-\frac{\theta B}{2}} V e^{\frac{\theta B}{2}}`.
