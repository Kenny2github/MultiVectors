Compute blades and multivectors in arbitrary-dimensional space.

Installation
------------

    pip install multivectors

Usage
-----

    >>> import math
    >>> from multivectors import x, y, z
    >>> v = 2*x + 3*y + 4*z
    >>> print(v.rotate(math.pi/2, x * y))
    (-3.00x + 2.00y + 4.00z)

For more see `the docs <https://multivectors.rtfd.io>`_