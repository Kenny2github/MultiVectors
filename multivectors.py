"""
# multivectors
Compute blades and multivectors in arbitrary-dimensional space.

## Installation
```
pip install multivectors
```

## Usage
```python
>>> import math
>>> from multivectors import x, y, z
>>> v = 2*x + 3*y + 4*z
>>> print(v.rotate(math.pi/2, x * y))
(-3.00x + 2.00y + 4.00z)
```

For more see [the docs](https://github.com/Kenny2github/MultiVectors/wiki)
"""
from __future__ import annotations
from itertools import combinations
import math
from numbers import Real
from typing import Dict, List, Tuple, Union

__all__ = [
    'Blade',
    'MultiVector',
    '_',
    'x',
    'y',
    'z',
    'w'
]

__version__ = '0.1.0'

NAMES = 'xyzw'

def merge(arr: List[int], left: List[int], right: List[int]) -> int:
    """Perform a merge of sorted lists and count the swaps."""
    i = j = count = 0
    ll = len(left)
    lr = len(right)
    while i < ll or j < lr:
        if i == ll:
            arr[i+j] = right[j]
            j += 1
        elif j == lr:
            arr[i+j] = left[i]
            i += 1
        elif left[i] <= right[j]:
            arr[i+j] = left[i]
            i += 1
        else:
            arr[i+j] = right[j]
            count += ll - i
            j += 1
    return count

def count_swaps(arr: List[int], copy: bool = True) -> int:
    """Count the number of swaps needed to sort a list.

    Examples:
    ```python
    >>> count_swaps([1, 3, 2, 5, 4])
    2
    >>> count_swaps([3, 2, 1])
    3

    ```
    """
    if copy:
        arr = list(arr)
    if len(arr) < 2:
        return 0
    m = (len(arr) + 1) // 2
    left = arr[:m]
    right = arr[m:]
    return (count_swaps(left, False)
            + count_swaps(right, False)
            + merge(arr, left, right))

def names_to_idxs(name: str, raise_on_invalid_chars: bool = False) -> List[int]:
    """Convert swizzled basis vector names into generalized basis indices.

    Examples:
    ```python
    >>> names_to_idxs('xyw')
    [0, 1, 3]
    >>> names_to_idxs('e_1e2_z')
    [0, 1, 2]
    >>> names_to_idxs('_')
    []

    ```
    """
    if name.startswith('__'):
        # fail on magic attributes immediately
        raise AttributeError
    idxs: List[int] = []
    i = 0
    nl = len(name)
    while i < nl:
        c = name[i]
        if c in NAMES:
            idxs.append(NAMES.index(c))
            i += 1
        elif c == 'e':
            i += 1
            j = i
            subname = []
            while j < nl and name[j] not in (NAMES + 'e') and not subname:
                if '0' <= name[j] <= '9':
                    subname.append(name[j])
                elif raise_on_invalid_chars:
                    raise AttributeError
                j += 1
            try:
                idxs.append(int(''.join(subname)) - 1)
            except ValueError:
                raise AttributeError('empty eN notation: %r' % name) from None
            i = j
        elif c == '_':
            i += 1
        elif raise_on_invalid_chars:
            raise AttributeError
        else:
            i += 1
    return idxs

def idxs_to_idxs(idxs: Index) -> List[int]:
    """Convert multiple possible ways to specify multiple indices.

    Examples:
    ```python
    >>> idxs_to_idxs(slice(None, 5, None))
    [0, 1, 2, 3, 4]
    >>> idxs_to_idxs((1, 3, 4))
    [1, 3, 4]
    >>> idxs_to_idxs(1)
    [1]

    ```
    """
    if isinstance(idxs, int):
        return [idxs]
    if isinstance(idxs, slice):
        if idxs.stop is None:
            raise TypeError('cannot have infinite bases')
        return list(range(*idxs.indices(idxs.stop)))
    return list(idxs)

class _BladeGetattr(type):
    def __getattr__(self, name: str) -> Blade:
        return Blade(*names_to_idxs(name))

    def __getitem__(self, idxs: Index) -> Blade:
        return Blade(*idxs_to_idxs(idxs))

class Blade(metaclass=_BladeGetattr):
    """Zero or more basis vectors multiplied, with a magnitude.

    *bases: 0-indices of individual basis vectors. Blade(0, 1) is the
        pair of x-hat * y-hat, aka e1 * e2
    scalar: Real number to scale the blade by.

    With zero bases, represents a scalar. In this capacity, functions
        pretty much entirely as a float but promotes floats to 0-blades.
    With one basis, represents a (scaled) basis vector, such as x-hat.
    With two or more bases, represents a k-blade, such as x-hat*y-hat.

    Basis vector names can be swizzled on the class itself:
    ```python
    >>> Blade.xy == Blade(0, 1)
    True
    >>> Blade.e1_e3 == Blade(0, 2)
    True
    >>> Blade._ == Blade()
    True

    ```
    And indices can be combined:
    ```python
    >>> Blade[1, 3] == Blade(1, 3)
    True
    >>> Blade[1:3] == Blade(1, 2)
    True

    ```
    """
    bases: Tuple[int, ...]
    scalar: Real

    _insts = {}

    @property
    def grade(self) -> int:
        """The k of "k-blade".

        ```python
        >>> (3 * Blade._).grade # scalars
        0
        >>> (2 * Blade.x).grade # vectors
        1
        >>> (4 * Blade.yz).grade # bivectors
        2

        ```
        """
        return len(self.bases)

    @property
    def terms(self) -> Tuple[Blade]:
        """For compatibility with MultiVector.

        ```python
        >>> Blade.x.terms
        (1.0 * Blade.x,)

        ```
        """
        return (self,)

    def __new__(cls, *bases: int, scalar: Real = 1.0) -> Simple:
        key = cls.condense_bases(bases, scalar)
        #if len(key[0]) == 0:
        #    return key[1]
        if key not in cls._insts:
            cls._insts[key] = super().__new__(cls)
        return cls._insts[key]

    def __init__(self, *bases: int, scalar: Real = 1.0) -> None:
        self.bases, self.scalar = self.condense_bases(bases, scalar)

    @classmethod
    def check(cls, obj: Simple) -> Blade:
        """Check if an object is a primitive real number.
        If so, promote it to a 0-blade. Otherwise, just return the object.
        This is mostly used internally.

        ```python
        >>> isinstance(2, Blade)
        False
        >>> isinstance(Blade.x, Blade)
        True
        >>> isinstance(Blade.check(2), Blade)
        True
        >>> isinstance(Blade.check(Blade.x), Blade)
        True

        ```
        """
        if isinstance(obj, Real):
            return cls(scalar=obj)
        return obj

    @staticmethod
    def condense_bases(bases: Tuple[int, ...], scalar: Real = 1.0) \
            -> Tuple[Tuple[int, ...], Real]:
        """Normalize a sequence of bases, modifying the scalar as necessary.

        bases: The tuple of basis indices.
        scalar: Real number that will scale the resulting blade.

        Returns: a 2-tuple of normalized bases and the modified scalar.

        Examples:
        ```python
        >>> Blade.condense_bases((1, 1, 2, 1, 2), 2.0)
        ((1,), -2.0)
        >>> Blade.condense_bases((1, 2, 1, 2), 1.5)
        ((), -1.5)
        >>> Blade.condense_bases((2, 1, 3, 2, 3, 3), 1.0)
        ((1, 3), 1.0)

        ```
        """
        bases = list(bases)
        if count_swaps(bases) % 2 == 1:
            scalar *= -1
        bases.sort()
        bases = [basis for basis in set(bases)
                 for _ in range(bases.count(basis) % 2)]
        return tuple(bases), scalar

    def __repr__(self) -> str:
        """Return the representation of this blade.

        If Blade is in the global namespace, this produces
        eval()-able code that gets the original object back.

        ```python
        >>> repr(Blade.xyzw)
        '1.0 * Blade.xyzw'
        >>> repr(Blade.e1e2e3e4e5)
        '1.0 * Blade(0, 1, 2, 3, 4)'

        ```
        """
        if not self.bases:
            return repr(self.scalar)
        if max(self.bases) > 3: # i=3 => e4 = w
            r = 'Blade(%s)' % ', '.join(map(repr, self.bases))
        else:
            r = 'Blade.' + ''.join(NAMES[i] for i in self.bases)
        return '%r * %s' % (self.scalar, r)

    def __str__(self) -> str:
        """Return the pretty representation of this blade.

        This prefers basis names over eval()-ability.

        ```python
        >>> str(Blade.xyzw)
        '1.00xyzw'
        >>> str(Blade.e1e2e3e4e5)
        '1.00e1e2e3e4e5'

        ```
        """
        if not self.bases:
            return str(self.scalar)
        if max(self.bases) > 3:
            r = ''.join(f'e{i+1}' for i in self.bases)
        else:
            r = ''.join(NAMES[i] for i in self.bases)
        return '%.2f%s' % (self.scalar, r)

    # Binary operators

    def __add__(self, other: MVV) -> MV:
        """Add a blade and another object.

        Returns: 0-blade if this blade has grade 0 and the object is scalar.
        Returns: k-blade with scalars summed if the object is a k-blade.
        Returns: MultiVector summing this object and the other otherwise.

        ```python
        >>> Blade._ * 2 + 1
        3.0
        >>> isinstance(Blade._ * 2 + 1, Blade)
        True
        >>> 2 * Blade.xy + 3 * Blade.xy
        5.0 * Blade.xy
        >>> Blade.x + 2 * Blade.y
        (1.0 * Blade.x + 2.0 * Blade.y)
        >>> Blade.x + Blade.xy
        (1.0 * Blade.x + 1.0 * Blade.xy)

        ```
        """
        if isinstance(other, Real) and self.bases == ():
            return Blade(scalar=self.scalar + other)
        if isinstance(other, Blade) and self.bases == other.bases:
            return Blade(*self.bases, scalar=self.scalar + other.scalar)
        if isinstance(other, _MVV):
            return MultiVector.from_terms(self, other)
        return NotImplemented

    def __radd__(self, other: Real) -> Blade:
        """Support adding blades on the right side of scalars.

        ```python
        >>> 1 + 2 * Blade.x
        (1.0 + 2.0 * Blade.x)
        >>> 2 + 3 * Blade._
        5.0
        >>> isinstance(2 + 3 * Blade._, Blade)
        True

        """
        if isinstance(other, Real):
            return self + other # scalar addition is commutative
        return NotImplemented

    def __sub__(self, other: MVV) -> MV:
        """Subtracting is adding the negation.

        ```python
        >>> Blade._ * 2 - 1
        1.0
        >>> isinstance(Blade._ * 2 - 1, Blade)
        True
        >>> 2 * Blade.xy - 3 * Blade.xy
        -1.0 * Blade.xy
        >>> Blade.x - 2 * Blade.y
        (1.0 * Blade.x + -2.0 * Blade.y)
        >>> Blade.x - Blade.xy
        (1.0 * Blade.x + -1.0 * Blade.xy)

        ```
        """
        return self + (-other)

    def __rsub__(self, other: Real) -> Blade:
        """Support subtracting blades from scalars.

        ```python
        >>> 1 - 2 * Blade.x
        (1.0 + -2.0 * Blade.x)
        >>> 2 - 3 * Blade._
        -1.0
        >>> isinstance(2 - 3 * Blade._, Blade)
        True

        ```
        """
        return other + (-self)

    def __mul__(self, other: Simple) -> Blade:
        """Multiply a blade and another object.

        Returns: k-blade scaled by object if object is scalar.
        Returns: The geometric product, if both objects are blades.

        ```python
        >>> Blade.x * 3 * 2
        6.0 * Blade.x
        >>> Blade.x * Blade.y
        1.0 * Blade.xy
        >>> Blade.z * 2 * Blade.x
        -2.0 * Blade.xz

        ```
        """
        if isinstance(other, Real):
            return Blade(*self.bases, scalar=self.scalar * other)
        if isinstance(other, Blade):
            return Blade(*self.bases, *other.bases,
                         scalar=self.scalar * other.scalar)
        return NotImplemented

    def __rmul__(self, other: Real) -> Blade:
        """Support multiplying blades on the right side of scalars.

        ```python
        >>> 2 * Blade.x * 3
        6.0 * Blade.x

        ```
        """
        if isinstance(other, Real):
            return self * other # scalar multiplication is commutative
        return NotImplemented

    def __matmul__(self, other: Simple) -> Simple:
        """Get the inner (dot) product of two objects.

        Returns: The scaled blade when dotting with a scalar.
        Returns: The inner product of the two blades, which can be scalar 0.

        ```python
        >>> Blade.x @ 3
        3.0 * Blade.x
        >>> Blade.x @ Blade.y
        0.0
        >>> Blade.x @ Blade.x
        1.0

        ```
        """
        if not isinstance(other, _Simple):
            return NotImplemented
        return (self * other) % abs(self.grade - Blade.check(other).grade)

    dot = inner = __matmul__

    def __rmatmul__(self, other: Real) -> Blade:
        """Support dotting blades on the right side of scalars.

        ```python
        >>> 2 @ Blade.x
        2.0 * Blade.x

        ```
        """
        if not isinstance(other, Real):
            return NotImplemented
        return (other * self) % abs(Blade.check(other).grade - self.grade)

    def __truediv__(self, other: Simple) -> Blade:
        """Divide two objects.

        Returns: The scaled blade when dividing by a scalar.
        Returns: The geometric product of this blade and the other's inverse.

        ```python
        >>> Blade.x / 2
        0.5 * Blade.x
        >>> # dividing basis vectors is the same as multiplying them
        >>> # since 1/x = x/||x||^2 = x
        >>> Blade.x / Blade.y
        1.0 * Blade.xy

        ```
        """
        if Blade.check(other).bases == ():
            return self * (1 / float(other))
        if isinstance(other, Blade):
            return self * other / (other * other)
        return NotImplemented

    def __rtruediv__(self, other: Real) -> Blade:
        """Divide a scalar by a blade.

        ```python
        >>> 1 / Blade.x
        1.0 * Blade.x

        ```
        """
        return other * self / (self * self)

    def __mod__(self, other: int) -> Simple:
        """The choose operator - returns the sum of all blades of grade k.
        When applied to a blade, returns the blade if k is the blade's grade.
        Returns scalar 0 otherwise.

        ```python
        >>> (2 * Blade.x) % 1
        2.0 * Blade.x
        >>> (2 * Blade.x) % 2
        0.0
        >>> (2 * Blade.xy) % 2
        2.0 * Blade.xy
        >>> Blade._ % 0
        1.0

        ```
        """
        return self if self.grade == other else 0.0

    choose = __mod__

    def __pow__(self, other: int, modulo: int = None) -> Simple:
        """A blade raised to an integer power.
        A ** n = A * A * A * ... * A n times.
        A ** -n = 1 / (A ** n)
        The optional ternary power uses the choose operator afterwards.

        ```python
        >>> Blade.xy ** 3
        -1.0 * Blade.xy
        >>> Blade.y ** -5
        1.0 * Blade.y
        >>> pow(Blade.xyz, 3, 3)
        -1.0 * Blade.xyz

        ```
        """
        if not isinstance(other, int):
            if not self.bases:
                return float(self) ** other
            return NotImplemented
        result = 1
        for _ in range(abs(other)):
            result *= self
        if other < 0:
            result = 1 / result
        if modulo is not None:
            if not isinstance(modulo, int):
                return NotImplemented
            return result % modulo
        return result

    def __rpow__(self, other: Real) -> Simple:
        """A real number raised to a blade power.
        x ** A = e ** (A ln x) = e ** (I * a ln x)
        = cos(a ln x) + sin(a ln x) * I

        ```python
        >>> from math import e, pi
        >>> from cmath import isclose
        >>> round(e ** (pi / 4 * Blade.xy), 2)
        (0.71 + 0.71 * Blade.xy)
        >>> # in 2D, xy is isomorphic to i
        >>> isclose(e ** (pi * Blade.xy), -1)
        True
        >>> isclose(e ** (pi * 1j), -1)
        True

        ```
        """
        if not isinstance(other, Real):
            return NotImplemented
        # Assume this vector is i
        theta = self.scalar * math.log(other)
        return math.cos(theta) + Blade(*self.bases, scalar=math.sin(theta))

    def __xor__(self, other: Simple) -> Simple:
        """Get the outer (wedge) product of two objects.
        WARNING: Operator precedence puts ^ after +!
        Make sure to put outer products in parentheses, like this:
        ``u * v == u @ v + (u ^ v)``

        Returns: The scaled blade when wedging with a scalar.
        Returns: The outer product of the two blades.

        ```python
        >>> Blade.x ^ 2
        2.0 * Blade.x
        >>> Blade.x ^ Blade.y
        1.0 * Blade.xy
        >>> 2 * Blade.x ^ Blade.y ^ 3 * Blade.w
        6.0 * Blade.xyw

        ```
        """
        if not isinstance(other, _Simple):
            return NotImplemented
        return (self * other) % (self.grade + Blade.check(other).grade)

    wedge = outer = __xor__

    def __rxor__(self, other: Real) -> Blade:
        """Support wedging blades on the right side of scalars.

        ```python
        >>> 2 ^ Blade.x
        2.0 * Blade.x

        ```
        """
        if not isinstance(other, Real):
            return NotImplemented
        return (other * self) % (Blade.check(other).grade + self.grade)

    # Unary operators

    def __neg__(self) -> Blade:
        """The negation of a blade is the blade with its scalar negated.

        ```python
        >>> -Blade.x
        -1.0 * Blade.x

        ```
        """
        return self * -1

    def __abs__(self) -> Union[Real, complex]:
        """Get the magnitude of a blade.

        ```python
        >>> from cmath import isclose
        >>> abs(Blade.x)
        1.0
        >>> isclose(abs(2 * Blade.x * 3 * Blade.y), 6j)
        True

        ```
        """
        # choosing the +ve square root because absolute value is +ve
        return (self * self) ** .5

    magnitude = __abs__

    def __invert__(self) -> Blade:
        """A normalized blade is one with a magnitude of 1.

        ```python
        >>> ~Blade.x # no effect on bases
        1.0 * Blade.x
        >>> ~(3 * Blade.xy)
        1.0 * Blade.xy

        ```
        """
        return Blade(*self.bases)

    normalize = __invert__

    def __complex__(self) -> complex:
        """Convert a scalar blade to a complex number.

        ```python
        >>> complex(2 * Blade._)
        (2+0j)
        >>> complex(3 * Blade.x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * Blade.x) to complex

        ```
        """
        if self.bases == ():
            return complex(self.scalar)
        raise TypeError(
            'cannot convert non-scalar blade (%r) to complex' % self)

    def __int__(self) -> int:
        """Convert a scalar blade to an integer.

        ```python
        >>> int(2 * Blade._)
        2
        >>> int(3 * Blade.x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * Blade.x) to int

        ```
        """
        if self.bases == ():
            return int(self.scalar)
        raise TypeError('cannot convert non-scalar blade (%r) to int' % self)

    def __float__(self) -> float:
        """Convert a scalar blade to a float.
        NOTE: This converts scalar blades only! To get the magnitude
        of a blade, use ``abs(blade)``.

        ```python
        >>> float(2 * Blade._)
        2.0
        >>> float(3 * Blade.x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * Blade.x) to float

        ```
        """
        if self.bases == ():
            return float(self.scalar)
        raise TypeError('cannot convert non-scalar blade (%r) to float' % self)

    def __index__(self) -> int:
        """Convert a scalar blade to an index.

        ```python
        >>> 'abcde'[2 * Blade._]
        'c'
        >>> 'abcde'[3 * Blade.x]
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * Blade.x) to int
        >>> 'abcde'[1.5 * Blade._]
        Traceback (most recent call last):
            ...
        TypeError: cannot losslessly convert float 1.5 to index

        ```
        """
        if int(self) == float(self):
            return int(self)
        raise TypeError('cannot losslessly convert float %r to index' % self)

    def __round__(self, ndigits: int = None) -> Blade:
        """Round the scalar value of a blade.

        ```python
        >>> round(1.26 * Blade.xy, 1)
        1.3 * Blade.xy
        >>> round(3.7 * Blade.zw)
        4 * Blade.zw

        ```
        """
        return Blade(*self.bases, scalar=round(self.scalar, ndigits))

    def __trunc__(self) -> Blade:
        """Truncate the scalar value of a blade.

        ```python
        >>> import math
        >>> math.trunc(1.2 * Blade.xy)
        1 * Blade.xy
        >>> math.trunc(-1.2 * Blade.zw)
        -1 * Blade.zw

        ```
        """
        return Blade(*self.bases, scalar=math.trunc(self.scalar))

    def __floor__(self) -> Blade:
        """Floor the scalar value of a blade.

        ```python
        >>> import math
        >>> math.floor(1.2 * Blade.xy)
        1 * Blade.xy
        >>> math.floor(-1.2 * Blade.zw)
        -2 * Blade.zw

        ```
        """
        return Blade(*self.bases, scalar=math.floor(self.scalar))

    def __ceil__(self) -> Blade:
        """Ceiling the scalar value of a blade.

        ```python
        >>> import math
        >>> math.ceil(1.2 * Blade.xy)
        2 * Blade.xy
        >>> math.ceil(-1.2 * Blade.zw)
        -1 * Blade.zw

        ```
        """
        return Blade(*self.bases, scalar=math.ceil(self.scalar))

class MultiVector:
    """Two or more incompatible blades, summed together.

    The bare constructor is not meant for regular use.
    Use the factory ``MultiVector.from_terms()`` instead.

    Basis vector names can be swizzled on instances:
    ```python
    >>> (Blade.x + Blade.y).x
    1.0
    >>> (Blade.x * Blade.y + Blade.z).xy
    1.0
    >>> (Blade.x + Blade.y).e3
    0.0

    ```
    And indices can be combined:
    ```python
    >>> (Blade.x + Blade.y)[0]
    1.0
    >>> (Blade.x * Blade.y + Blade.z)[0, 1]
    1.0
    >>> (Blade.x + Blade.y)[2]
    0.0

    ```
    """

    termdict: Dict[Tuple[int, ...], Real]

    @property
    def terms(self) -> Tuple[Blade, ...]:
        """Get a sequence of blades comprising this multivector.

        ```python
        >>> (Blade.x + Blade.y).terms
        (1.0 * Blade.x, 1.0 * Blade.y)
        >>> ((Blade.x + Blade.y) * (Blade.z + Blade.w)).terms
        (1.0 * Blade.xz, 1.0 * Blade.xw, 1.0 * Blade.yz, 1.0 * Blade.yw)

        ```
        """
        keys = sorted(self.termdict.keys(), key=lambda b: (len(b), b))
        return tuple(Blade(*key, scalar=self.termdict[key]) for key in keys)

    def __init__(self, termdict: Dict[Tuple[int, ...], Real]):
        self.termdict = termdict.copy()

    @classmethod
    def from_terms(cls, *terms: MV) -> MV:
        """Create a multivector from a sequence of terms.
        This may return something other than a MultiVector
        if it is the only term or there are no terms.

        ```python
        >>> MultiVector.from_terms(Blade.x, Blade.y)
        (1.0 * Blade.x + 1.0 * Blade.y)
        >>> MultiVector.from_terms()
        0.0
        >>> MultiVector.from_terms(Blade.z)
        1.0 * Blade.z
        >>> MultiVector.from_terms(2 * Blade.x, Blade.x)
        3.0 * Blade.x

        ```
        """
        terms: List[Simple] = [
            t for term in terms for t in
            (term.terms if isinstance(term, MultiVector) else (term,))
        ]
        termdict: Dict[Tuple[int, ...], List[Real]] = {}
        for t in terms:
            if abs(t) == 0:
                continue
            if isinstance(t, Real):
                t = Blade(scalar=t)
            termdict.setdefault(t.bases, []).append(t.scalar)
        d = {b: math.fsum(s) for b, s in termdict.items()}
        for key in [k for k in d if d[k] == 0]:
            del d[key]
        inst = cls(d)
        ts = list(inst.terms)
        if not ts:
            return 0.0
        if len(ts) == 1:
            return ts[0]
        return inst

    def __getattr__(self, name: str) -> Union[Real, Tuple[Real, ...]]:
        """Support basis name swizzling."""
        return self.termdict.get(tuple(names_to_idxs(name)), 0.0)

    def __getitem__(self, idxs: Index) -> Union[Real, Tuple[Real]]:
        """Support index swizzling."""
        return self.termdict.get(tuple(idxs_to_idxs(idxs)), 0.0)

    def __repr__(self) -> str:
        """Return a representation of this multivector.
        Depending on the global namespace, this may be eval()-able.

        ```python
        >>> repr(Blade.x + Blade.y)
        '(1.0 * Blade.x + 1.0 * Blade.y)'
        >>> repr(Blade.yz - Blade.xw)
        '(-1.0 * Blade.xw + 1.0 * Blade.yz)'

        ```
        """
        return '(' + ' + '.join(map(repr, self.terms)) + ')'

    def __str__(self) -> str:
        """Return a representation of this multivector suited for showing.

        ```python
        >>> str(Blade.x + Blade.y)
        '(1.00x + 1.00y)'
        >>> str(Blade.yz - Blade.xw)
        '(-1.00xw + 1.00yz)'

        ```
        """
        return '(' + ' + '.join(map(str, self.terms)) + ')'

    # Binary operators

    def __add__(self, other: MVV) -> MV:
        """Add a multivector and another object.

        Returns: The sum of the terms of this multivector and the other.
        Returns: The other object added to this multivector's terms.

        ```python
        >>> (Blade.x + Blade.z) + (Blade.y + Blade.w)
        (1.0 * Blade.x + 1.0 * Blade.y + 1.0 * Blade.z + 1.0 * Blade.w)
        >>> (Blade.x + Blade.z) + Blade.y
        (1.0 * Blade.x + 1.0 * Blade.y + 1.0 * Blade.z)

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*self.terms, *other.terms)
        if isinstance(other, _Simple):
            return self.from_terms(*self.terms, other)
        return NotImplemented

    def __radd__(self, other: MVV) -> MV:
        """Support adding multivectors on the right side of objects.

        ```python
        >>> 1 + (Blade.x + Blade.z)
        (1.0 + 1.0 * Blade.x + 1.0 * Blade.z)
        >>> Blade.x + (Blade.y + Blade.z)
        (1.0 * Blade.x + 1.0 * Blade.y + 1.0 * Blade.z)

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*other.terms, *self.terms)
        if isinstance(other, _Simple):
            return self.from_terms(other, *self.terms)
        return NotImplemented

    def __sub__(self, other: MVV) -> MV:
        """Subtracting is adding the negation.

        ```python
        >>> (Blade.x + Blade.z) - (Blade.y + Blade.w)
        (1.0 * Blade.x + -1.0 * Blade.y + 1.0 * Blade.z + -1.0 * Blade.w)
        >>> (Blade.x + Blade.z) - Blade.y
        (1.0 * Blade.x + -1.0 * Blade.y + 1.0 * Blade.z)

        ```
        """
        return self + (-other)

    def __rsub__(self, other: MVV) -> MV:
        """Support subtracting multivectors from objects.

        ```python
        >>> 1 - (Blade.x + Blade.z)
        (1.0 + -1.0 * Blade.x + -1.0 * Blade.z)
        >>> Blade.x - (Blade.y + Blade.z)
        (1.0 * Blade.x + -1.0 * Blade.y + -1.0 * Blade.z)

        ```
        """
        return other + (-self)

    def __mul__(self, other: MVV) -> MV:
        """Multiply a multivector and another object.

        Returns: (a+b)*(c+d)=a*c+a*d+b*c+b*d for multivectors (a+b) and (c+d)
        Returns: (a+b)*v = a*v + b*v for multivector (a+b) and simple v

        ```python
        >>> (Blade.x + Blade.y) * (Blade.z + Blade.w)
        (1.0 * Blade.xz + 1.0 * Blade.xw + 1.0 * Blade.yz + 1.0 * Blade.yw)
        >>> (Blade.x + Blade.y) * 3
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> (Blade.x + Blade.y) * Blade.x
        (1.0 + -1.0 * Blade.xy)

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*(a * b for b in other.terms
                                     for a in self.terms))
        if isinstance(other, _Simple):
            return self.from_terms(*(t * other for t in self.terms))
        return NotImplemented

    def __rmul__(self, other: Simple) -> MV:
        """Support multiplying multivectors on the right side of simples.

        ```python
        >>> 3 * (Blade.x + Blade.y)
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> Blade.y * (Blade.x + Blade.y)
        (1.0 + -1.0 * Blade.xy)

        ```
        """
        if isinstance(other, _Simple):
            return self.from_terms(*(other * t for t in self.terms))
        return NotImplemented

    def __matmul__(self, other: MVV) -> MV:
        """Get the inner (dot) product of two objects.

        Returns: (a+b)@(c+d)=a@c+a@d+b@c+b@d for multivectors (a+b) and (c+d)
        Returns: (a+b)@v = a@v + b@v for multivector (a+b) and simple v

        ```python
        >>> (2 * Blade.x + 3 * Blade.y) @ (4 * Blade.x + 5 * Blade.y)
        23.0
        >>> (2 * Blade.xy) @ (3 * Blade.yz)
        0.0
        >>> (Blade.x + Blade.y) @ 3
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> (Blade.x + Blade.y) @ Blade.x
        1.0

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*(
                Blade.check(a) @ Blade.check(b)
                for b in other.terms for a in self.terms))
        if isinstance(other, _Simple):
            return self.from_terms(*(
                Blade.check(t) @ Blade.check(other)
                for t in self.terms))
        return NotImplemented

    dot = inner = __matmul__

    def __rmatmul__(self, other: Simple) -> MV:
        """Support dotting multivectors on the right hand side.

        Returns: v@(a+b) = v@a + v@b for multivector (a+b) and simple v

        ```python
        >>> 3 @ (Blade.x + Blade.y)
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> Blade.x @ (Blade.x + Blade.y)
        1.0

        ```
        """
        if isinstance(other, _Simple):
            return self.from_terms(*(
                Blade.check(other) @ Blade.check(t)
                for t in self.terms))
        return NotImplemented

    def __truediv__(self, other: Simple) -> MV:
        """Divide two objects.

        Returns: (a+b)/v = a/v + b/v for multivector (a+b) and simple v

        ```python
        >>> (6 * Blade.x + 9 * Blade.y) / 3
        (2.0 * Blade.x + 3.0 * Blade.y)
        >>> (6 * Blade.x + 9 * Blade.y) / (3 * Blade.x)
        (2.0 + -3.0 * Blade.xy)

        ```
        """
        if isinstance(other, _Simple):
            return self.from_terms(*(t / other for t in self.terms))
        return NotImplemented

    def __mod__(self, grade: int) -> MVV:
        """The choose operator - returns the sum of all blades of grade k.

        Examples:
        ```python
        >>> (1 + 2 * Blade.x + 3 * Blade.y + 4 * Blade.xy) % 1
        (2.0 * Blade.x + 3.0 * Blade.y)
        >>> (1 + 2 + 3 * Blade.x + 4 * Blade.xy + 5 * Blade.yz) % 2
        (4.0 * Blade.xy + 5.0 * Blade.yz)
        >>> (1 + 2 + 3 * Blade.x + 4 * Blade.xy + 5 * Blade.yz) % 0
        3.0
        """
        if not isinstance(grade, int):
            return NotImplemented
        bs = set(basis for bases in self.termdict.keys() for basis in bases)
        return sum(Blade(*bases, scalar=self[bases])
                   for bases in combinations(bs, grade))

    choose = __mod__

    def __pow__(self, other: int, modulo: int = None) -> MVV:
        """A multivector raised to an integer power.
        V ** n = V * V * V * ... * V n times.
        The optional ternary power uses the choose operator afterwards.
        Note that unlike blades, multivectors do not have inverses, so
        negative powers are not supported.

        ```python
        >>> (Blade.x + Blade.y) ** 3
        (2.0 * Blade.x + 2.0 * Blade.y)
        >>> pow(Blade.x + Blade.yz, 3, 2)
        2.0 * Blade.yz

        """
        if not isinstance(other, int) or other < 0:
            return NotImplemented
        result = 1
        for _ in range(other):
            result *= self
        if modulo is not None:
            if not isinstance(modulo, int):
                return NotImplemented
            return result % modulo
        return result

    def __rpow__(self, other: Real) -> MVV:
        """A real number raised to a multivector power.
        x ** V = Product of x ** Vi = Product of e ** (Vi ln x)
        = Product of (cos(ai ln x) + sin(ai ln x) * Ii)

        ```python
        >>> from math import e, pi
        >>> round(e ** (pi/4 * (Blade.x + Blade.y)), 2)
        (0.5 + 0.5 * Blade.x + 0.5 * Blade.y + 0.5 * Blade.xy)
        >>> round(e ** (pi * (Blade.x + Blade.y)), 2)
        1.0

        ```
        """
        if not isinstance(other, Real):
            return NotImplemented
        result = 1
        for t in self.terms:
            result *= other ** t
        return result

    def __xor__(self, other: MVV) -> MV:
        """Get the outer (wedge) product of two objects.
        WARNING: Operator precedence puts ^ after +!
        Make sure to put outer products in parentheses, like this:
        ``u * v == u @ v + (u ^ v)``

        Returns: (a+b)^(c+d)=(a^c)+(a^d)+(b^c)+(b^d) for (a+b), (c+d)
        Returns: (a+b)^v = (a^v) + (b^v) for multivector (a+b) and simple v

        ```python
        >>> (2 * Blade.x + 3 * Blade.y) ^ (4 * Blade.x + 5 * Blade.y)
        -2.0 * Blade.xy
        >>> (2 * Blade.xy) ^ (3 * Blade.yz)
        0.0
        >>> (Blade.x + Blade.y) ^ 3
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> (Blade.x + Blade.y) ^ Blade.x
        -1.0 * Blade.xy

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*(
                Blade.check(a) ^ Blade.check(b)
                for b in other.terms for a in self.terms))
        if isinstance(other, _Simple):
            return self.from_terms(*(
                Blade.check(t) ^ Blade.check(other)
                for t in self.terms))
        return NotImplemented

    wedge = outer = __xor__

    def __rxor__(self, other: Simple) -> MV:
        """Support wedging multivectors on the right hand side.

        Returns: v@(a+b) = v@a + v@b for multivector (a+b) and simple v

        ```python
        >>> 3 ^ (Blade.x + Blade.y)
        (3.0 * Blade.x + 3.0 * Blade.y)
        >>> Blade.x ^ (Blade.x + Blade.y)
        1.0 * Blade.xy

        ```
        """
        if isinstance(other, _Simple):
            return self.from_terms(*(
                Blade.check(other) ^ Blade.check(t)
                for t in self.terms))
        return NotImplemented

    # Unary operators

    def __neg__(self) -> MultiVector:
        """The negation of a multivector is the negation of all its terms.

        ```python
        >>> -(2 * Blade.x + 3 * Blade.yz)
        (-2.0 * Blade.x + -3.0 * Blade.yz)

        ```
        """
        return self.from_terms(*(-t for t in self.terms))

    def __abs__(self) -> Real:
        """The magnitude of a multivector is the square root of
        the sum of the squares of its terms.

        ```python
        >>> abs(Blade.x + Blade.y + Blade.z + Blade.w)
        2.0
        >>> abs(2 ** 1.5 * Blade.x + 2 ** 1.5 * Blade.y)
        4.0

        ```
        """
        return math.fsum(t * t for t in self.terms) ** .5

    magnitude = __abs__

    def __invert__(self) -> MultiVector:
        """A normalized multivector is one scaled down by its magnitude.

        ```python
        >>> ~(Blade.x + Blade.y + Blade.z + Blade.w)
        (0.5 * Blade.x + 0.5 * Blade.y + 0.5 * Blade.z + 0.5 * Blade.w)
        >>> round(~(Blade.x + Blade.y), 3)
        (0.707 * Blade.x + 0.707 * Blade.y)

        ```
        """
        return self / abs(self)

    normalize = __invert__

    def __round__(self, ndigits: int = None) -> MultiVector:
        """Round the scalars of each component term of a multivector.

        ```python
        >>> round(1.7 * Blade.x + 1.2 * Blade.y)
        (2.0 * Blade.x + 1.0 * Blade.y)
        >>> round(0.15 * Blade.x + 0.05 * Blade.y, 1)
        (0.1 * Blade.x + 0.1 * Blade.y)

        ```
        """
        return MultiVector.from_terms(*(round(t, ndigits) for t in self.terms))

    def __trunc__(self) -> MultiVector:
        """Truncate the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> math.trunc(1.7 * Blade.x + 1.2 * Blade.y)
        (1.0 * Blade.x + 1.0 * Blade.y)
        >>> math.trunc(-1.7 * Blade.x - 1.2 * Blade.y)
        (-1.0 * Blade.x + -1.0 * Blade.y)

        ```
        """
        return MultiVector.from_terms(*(math.trunc(t) for t in self.terms))

    def __floor__(self) -> MultiVector:
        """Floor the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> math.floor(1.7 * Blade.x + 1.2 * Blade.y)
        (1.0 * Blade.x + 1.0 * Blade.y)
        >>> math.floor(-1.7 * Blade.x - 1.2 * Blade.y)
        (-2.0 * Blade.x + -2.0 * Blade.y)

        ```
        """
        return MultiVector.from_terms(*(math.floor(t) for t in self.terms))

    def __ceil__(self) -> MultiVector:
        """Ceiling the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> math.ceil(1.7 * Blade.x + 1.2 * Blade.y)
        (2.0 * Blade.x + 2.0 * Blade.y)
        >>> math.ceil(-1.7 * Blade.x - 1.2 * Blade.y)
        (-1.0 * Blade.x + -1.0 * Blade.y)

        ```
        """
        return MultiVector.from_terms(*(math.ceil(t) for t in self.terms))

    # Actual methods

    def rotate(self, angle: Real, plane: Blade) -> MultiVector:
        """Rotate this multivector by angle in rads around the blade plane.

        angle: Angle to rotate by, in radians.
        plane: Blade representing basis plane to rotate through.

        ```python
        >>> from math import radians
        >>> round((3 * Blade.x + 2 * Blade.y + 4 * Blade.z).rotate(
        ...     radians(90), Blade.xy), 2)
        (-2.0 * Blade.x + 3.0 * Blade.y + 4.0 * Blade.z)
        >>> round((3 * Blade.x + 2 * Blade.y + 4 * Blade.z + 5 * Blade.w)
        ...     .rotate(radians(90), Blade.xyz), 2)
        (3.0 * Blade.x + 2.0 * Blade.y + 4.0 * Blade.z + -5.0 * Blade.xyzw)

        ```
        """
        power = plane * angle / 2
        return (math.e ** (-power)) * self * (math.e ** power)

    def angle_to(self, other: MultiVector) -> Real:
        """Get the angle between this multivector and another.
        NOTE: Only works for multivectors of singular grade 1!

        ```python
        >>> from math import degrees
        >>> math.degrees((Blade.x + Blade.y).angle_to(Blade.x - Blade.y))
        90.0
        >>> round(math.degrees((Blade.x + Blade.y + Blade.z + Blade.w).angle_to(
        ...     Blade.x - Blade.y - Blade.z - Blade.w)), 2)
        120.0

        ```
        """
        return math.acos((self @ other) / (abs(self) * abs(other)))

Simple = Union[Real, Blade]
_Simple = Real, Blade
MV = Union[Blade, MultiVector]
MVV = Union[Real, Blade, MultiVector]
_MVV = Real, Blade, MultiVector
Index = Union[int, Tuple[int, ...], slice]

def __getattr__(name: str) -> Blade:
    """Support module level swizzling of basis vectors.
    Unlike Blade.__getattr__, this rejects invalid characters in the name.
    """
    return Blade(*names_to_idxs(name, True))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
