As Applied in This Library
==========================

All of the concepts in :ref:`geo-intro` are applied in this library.

.. default-domain:: py

.. currentmodule:: multivectors

Blades
------

* :data:`multivectors.x`, :data:`y`, :data:`z`, and :data:`w` represent the first four basis vectors. Bases past 4D use :math:`\hat e_n` syntax: ``multivectors.e5`` or ``multivectors.e_5`` represent the 5th basis vector.
* Basis names can be *swizzled* on the module:

    .. code-block:: python

        >>> from multivectors import x, y, z, xyz
        >>> x * y * z == xyz
        True

* :data:`multivectors._` is a :math:`0`-blade: a scalar, but with :class:`MultiVector` type. It is differentiated from normal scalars when :func:`repr()`-d by surrounding parentheses:

    .. code-block:: python

        >>> from multivectors import _
        >>> 2 * _
        (2.0)

* The rules of arithmetic with blades as described above apply:

    .. code-block:: python

        >>> from multivectors import x, y, z
        >>> x * 2 + x * 3
        (5.0 * x)
        >>> x * 5 * y
        (5.0 * x*y)
        >>> (x * y) * z == x * (y * z)
        True

* You can query the grade of a blade (only if it is an actual blade!):

    .. code-block:: python

        >>> from multivectors import xy, yzw
        >>> xy.grade
        2
        >>> (2 * yzw).grade
        3

MultiVector
-----------

* A multivector is a sum of zero or more blades. All scalars and blades are multivectors, but not all multivectors are scalars or blades.
* Basis names can be swizzled on class **instances** to get the coefficient of that basis:

    .. code-block:: python

        >>> from multivectors import x, y
        >>> (x + 2*y).x
        1.0


* Basis indices can also be used:

    .. code-block:: python

        >>> from multivectors import xy, yz
        >>> (xy + 2*yz) % (0, 1)
        1.0


* Choosing by grade is supported:

    .. code-block:: python

        >>> from multivectors import x, y, xy, yz
        >>> V = 1 + 2*x + 3*y + 4*xy + 5*yz
        >>> V[0]
        (1.0)
        >>> V[1]
        (2.0 * x + 3.0 * y)
        >>> V[2]
        (4.0 * x*y + 5.0 * y*z)


* Getting the grade of a multivector will only work if if the multivector is a blade. Otherwise, you get :obj:`None`:

    .. code-block:: python

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

* The rules of arithmetic with multivectors as described above apply:

    .. code-block:: python

        >>> from multivectors import x, y, z
        >>> x * (y + z)
        (1.0 * x*y + 1.0 * x*z)
        >>> (y + z) * (2 * x)
        (-2.0 * x*y + -2.0 * x*z)
        >>> (1*x + 2*y) * (3*x + 4*y)
        (11.0 + -2.0 * x*y)

* The extra products apply too:

    .. code-block:: python

        >>> from multivectors import x, y, z
        >>> 1*3 + 2*4
        11
        >>> (1*x + 2*y) @ (3*x + 4*y)
        (11.0)
        >>> (1*5 - 2*4, 1*6 - 3*4, 2*6 - 3*5)
        (-3, -6, -3)
        >>> (1*x + 2*y + 3*z) ^ (4*x + 5*y + 6*z)
        (-3.0 * x*y + -6.0 * x*z + -3.0 * y*z)

* A convenience method is provided to rotate multivectors:

    .. code-block:: python

        >>> from math import radians
        >>> from multivectors import x, y, z, xz
        >>> round((3*x + 2*y + 4*z).rotate(radians(90), xz), 2)
        (-4.0 * x + 2.0 * y + 3.0 * z)
