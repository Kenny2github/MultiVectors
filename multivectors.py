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

For more see [the docs](https://github.com/Kenny2github/MultiVectors/blob/main/docs.md)
"""
from __future__ import annotations
import math
from numbers import Real
from typing import Dict, Iterable, List, Tuple, Type, Union

__all__ = [
    'Blade',
    'MultiVector',
    '_',
    'x',
    'y',
    'z',
    'w'
]

__version__ = '0.1.1'

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

def idxs_to_names(idxs: Index, sep='') -> str:
    """Convert indices to a swizzled name combination.

    Examples:
    ```python
    >>> idxs_to_names(slice(None, 5, None))
    'e1e2e3e4e5'
    >>> idxs_to_names((0, 1, 3))
    'xyw'
    >>> idxs_to_names(2)
    'z'

    ```
    """
    idxs = idxs_to_idxs(idxs)
    if max(idxs) > 3:
        return sep.join(f'e{i+1}' for i in idxs)
    return sep.join(NAMES[i] for i in idxs)

def condense_bases(bases: Tuple[int, ...], scalar: Scalar_a = None) \
        -> Tuple[Tuple[int, ...], Scalar_a]:
    """Normalize a sequence of bases, modifying the scalar as necessary.

    bases: The tuple of basis indices.
    scalar: Real number that will scale the resulting bases.

    Returns: a 2-tuple of normalized bases and the modified scalar.

    Examples:
    ```python
    >>> condense_bases((1, 1, 2, 1, 2), 2.0)
    ((1,), -2.0)
    >>> condense_bases((1, 2, 1, 2), 1.5)
    ((), -1.5)
    >>> condense_bases((2, 1, 3, 2, 3, 3), 1.0)
    ((1, 3), 1.0)

    ```
    """
    bases = list(bases)
    if scalar is None:
        scalar = Scalar_f('1')
    if count_swaps(bases) % 2:
        scalar *= Scalar_f('-1')
    bases.sort()
    bases = [basis for basis in set(bases)
                for _ in range(bases.count(basis) % 2)]
    return tuple(bases), scalar

class MultiVector:
    """The sum of one or more blades.

    The bare constructor is not meant for regular use.
    Use the factory ``MultiVector.from_terms()`` instead.

    Basis vector names can be swizzled on instances:
    ```python
    >>> from multivectors import x, y, z
    >>> (x + y).x
    1.0
    >>> (x*y + z).xy
    1.0
    >>> (x + y).e3
    0.0

    ```
    And indices can be combined:
    ```python
    >>> (x + y) % 0
    1.0
    >>> (x*y + z) % (0, 1)
    1.0
    >>> (x + y) % 2
    0.0

    ```
    """

    termdict: TermDict

    @property
    def grade(self) -> Union[int, None]:
        """The grade of this blade.

        Returns: None if this multivector is not a blade (one term)
        Returns: the number of different bases this blade consists of

        ```python
        >>> from multivectors import x, y, z
        >>> (x + y).grade
        >>> (z * 2 + z).grade
        1
        >>> (x*y*z).grade
        3

        ```
        """
        if len(self.termdict) != 1:
            return None
        (bases,) = self.termdict.keys()
        return len(bases)

    @property
    def terms(self) -> Tuple[MultiVector, ...]:
        """Get a sequence of blades comprising this multivector.

        ```python
        >>> from multivectors import x, y, z, w
        >>> (x + y).terms
        ((1.0 * x), (1.0 * y))
        >>> ((x + y) * (z + w)).terms
        ((1.0 * x*z), (1.0 * x*w), (1.0 * y*z), (1.0 * y*w))

        ```
        """
        items = sorted(self.termdict.items(), key=lambda i: (len(i[0]), i[0]))
        return tuple(MultiVector({key: value}) for key, value in items)

    def __init__(self, termdict: TermDict):
        self.termdict = {b: Scalar_f(s) for b, s in termdict.items()
                         if s != Scalar_f()}
        self.termdict = self.termdict or {(): Scalar_f()}

    @classmethod
    def from_terms(cls, *terms: SOV) -> MultiVector:
        """Create a multivector from a sequence of terms.

        ```python
        >>> from multivectors import x, y, z
        >>> MultiVector.from_terms(x, y)
        (1.0 * x + 1.0 * y)
        >>> MultiVector.from_terms()
        (0.0)
        >>> MultiVector.from_terms(z)
        (1.0 * z)
        >>> MultiVector.from_terms(2 * x, x)
        (3.0 * x)

        ```
        """
        termseq = [
            t for term in terms for t in
            (term if isinstance(term, MultiVector)
             else MultiVector({(): term})).termdict.items()
        ]
        termdict: Dict[Tuple[int, ...], List[Scalar_a]] = {}
        for bases, scalar in termseq:
            termdict.setdefault(bases, []).append(scalar)
        d = {b: Scalar_f(sum(s)) for b, s in termdict.items()}
        d = {b: s for b, s in d.items() if s != Scalar_f()} # discard zeroes
        return cls(d)

    def __getattr__(self, name: str) -> Scalar_a:
        """Support basis name swizzling."""
        return self.termdict.get(tuple(names_to_idxs(name)), Scalar_f())

    def __getitem__(self, grades: Index) -> MultiVector:
        """The choose operator - returns the sum of all blades of grade k.

        Examples:
        ```python
        >>> from multivectors import x, y, z
        >>> (1 + 2*x + 3*y + 4*x*y)[1]
        (2.0 * x + 3.0 * y)
        >>> (1 + 2 + 3*x + 4*x*y + 5*y*z)[2]
        (4.0 * x*y + 5.0 * y*z)
        >>> (1 + 2 + 3*x + 4*x*y + 5*y*z)[0]
        (3.0)
        >>> (1 + 2 + 3*x + 4*x*y + 5*y*z).choose(0)
        (3.0)
        >>> (1 + 2*x + 3*x*y)[:2]
        (1.0 + 2.0 * x)

        ```
        """
        grades = set(idxs_to_idxs(grades))
        return MultiVector({
            bases: scalar for bases, scalar in self.termdict.items()
            if len(bases) in grades})

    choose = __getitem__

    def __repr__(self) -> str:
        """Return a representation of this multivector.
        Depending on the global namespace, this may be eval()-able.

        ```python
        >>> from multivectors import x, y, z, w
        >>> repr(x + y)
        '(1.0 * x + 1.0 * y)'
        >>> repr(y*z - x*w)
        '(-1.0 * x*w + 1.0 * y*z)'

        ```
        """
        if self.grade is not None:
            ((bases, scalar),) = self.termdict.items()
            if bases == ():
                return '(%r)' % scalar
            return '(%r * %s)' % (scalar, idxs_to_names(bases, '*'))
        return '(' + ' + '.join(repr(t).strip('()') for t in self.terms) + ')'

    def __str__(self) -> str:
        """Return a representation of this multivector suited for showing.

        ```python
        >>> from multivectors import x, y, z, w
        >>> str(x + y)
        '(1.00x + 1.00y)'
        >>> str(y*z - x*w)
        '(-1.00xw + 1.00yz)'

        ```
        """
        if self.grade is not None:
            ((bases, scalar),) = self.termdict.items()
            if bases == ():
                return str(scalar)
            return '%.2f%s' % (scalar, idxs_to_names(bases))
        return '(' + ' + '.join(map(str, self.terms)) + ')'

    # Relational operators

    def __eq__(self, other: SOV) -> bool:
        """Compare equality of two objects.

        Returns: True if all terms of this multivector are equal to the other.
        Returns: True if this multivector is scalar and equals the other.
        Returns: False for all other cases or types.

        ```python
        >>> from multivectors import x, y
        >>> x + y == y + x
        True
        >>> x + 2*y == 2*x + y
        False

        ```
        """
        if not isinstance(other, MultiVector):
            if self.grade != 0:
                return False # this MV is not a scalar
            return self._ == other
        return self.termdict == other.termdict

    def __ne__(self, other: SOV) -> bool:
        """Compare inequality of two objects.

        Returns: False if all terms of this multivector are equal to the other.
        Returns: False if this multivector is scalar and equals the other.
        Returns: True for all other cases or types.

        ```python
        >>> from multivectors import x, y
        >>> x + y != y + x
        False
        >>> x + 2*y != 2*x + y
        True

        ```
        """
        return not (self == other)

    def __lt__(self, other: Scalar_a) -> bool:
        """Compare this blade less than an object.

        Returns: True if this is a scalar blade less than the scalar.
        Returns: NotImplemented for all other types.

        ```python
        >>> from multivectors import _, x
        >>> _ * 1 < 2
        True
        >>> _ * 2 < 1
        False
        >>> x * 1 < 2
        Traceback (most recent call last):
            ...
        TypeError: '<' not supported between instances of \
'MultiVector' and 'int'

        ```
        """
        if not isinstance(other, Scalar_t) or self.grade != 0:
            return NotImplemented
        return self._ < other

    def __gt__(self, other: Scalar_a) -> bool:
        """Compare this blade greater than an object.

        Returns: True if this is a scalar blade greater than the scalar.
        Returns: NotImplemented for all other types.

        ```python
        >>> from multivectors import _, x
        >>> _ * 1 > 2
        False
        >>> _ * 2 > 1
        True
        >>> x * 1 > 2
        Traceback (most recent call last):
            ...
        TypeError: '>' not supported between instances of \
'MultiVector' and 'int'

        ```
        """
        if not isinstance(other, Scalar_t) or self.grade != 0:
            return NotImplemented
        return self._ > other

    def __le__(self, other: Scalar_a) -> bool:
        """Compare this blade less than or equal to an object.

        Returns: True if this is a scalar blade not greater than the scalar.
        Returns: NotImplemented for all other types.

        ```python
        >>> from multivectors import _, x
        >>> _ * 1 <= 2
        True
        >>> _ * 2 <= 2
        True
        >>> _ * 2 <= 1
        False
        >>> x * 1 <= 2
        Traceback (most recent call last):
            ...
        TypeError: '<=' not supported between instances of \
'MultiVector' and 'int'

        ```
        """
        if not isinstance(other, Scalar_t) or self.grade != 0:
            return NotImplemented
        return self._ <= other

    def __ge__(self, other: Scalar_a) -> bool:
        """Compare this blade greater than or equal to an object.

        Returns: True if this is a scalar blade not less than the scalar.
        Returns: NotImplemented for all other types.

        ```python
        >>> from multivectors import _, x
        >>> _ * 1 >= 2
        False
        >>> _ * 2 >= 2
        True
        >>> _ * 2 >= 1
        True
        >>> x * 1 >= 2
        Traceback (most recent call last):
            ...
        TypeError: '>=' not supported between instances of \
'MultiVector' and 'int'

        ```
        """
        if not isinstance(other, Scalar_t) or self.grade != 0:
            return NotImplemented
        return self._ >= other

    # Binary operators

    def __add__(self, other: SOV) -> MultiVector:
        """Add a multivector and another object.

        Returns: The sum of the terms of this multivector and the other.
        Returns: The scalar added to this multivector's terms.

        ```python
        >>> from multivectors import x, y, z, w
        >>> (x + z) + (y + w)
        (1.0 * x + 1.0 * y + 1.0 * z + 1.0 * w)
        >>> (x + z) + y
        (1.0 * x + 1.0 * y + 1.0 * z)

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*self.terms, *other.terms)
        if isinstance(other, Scalar_t):
            return self.from_terms(*self.terms, other)
        return NotImplemented

    def __radd__(self, other: Scalar_a) -> MultiVector:
        """Support adding multivectors on the right side of objects.

        ```python
        >>> from multivectors import x, y, z
        >>> 1 + (x + z)
        (1.0 + 1.0 * x + 1.0 * z)
        >>> x + (y + z)
        (1.0 * x + 1.0 * y + 1.0 * z)

        ```
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*other.terms, *self.terms)
        if isinstance(other, Scalar_t):
            return self.from_terms(other, *self.terms)
        return NotImplemented

    def __sub__(self, other: SOV) -> MultiVector:
        """Subtracting is adding the negation.

        ```python
        >>> from multivectors import x, y, z, w
        >>> (x + z) - (y + w)
        (1.0 * x + -1.0 * y + 1.0 * z + -1.0 * w)
        >>> (x + z) - y
        (1.0 * x + -1.0 * y + 1.0 * z)

        ```
        """
        return self + (-other)

    def __rsub__(self, other: Scalar_a) -> MultiVector:
        """Support subtracting multivectors from objects.

        ```python
        >>> from multivectors import x, y, z, w
        >>> 1 - (x + z)
        (1.0 + -1.0 * x + -1.0 * z)
        >>> x - (y + z)
        (1.0 * x + -1.0 * y + -1.0 * z)

        ```
        """
        return other + (-self)

    def __mul__(self, other: SOV) -> MultiVector:
        """Multiply a multivector and another object.

        Returns: (a+b)*(c+d)=a*c+a*d+b*c+b*d for multivectors (a+b) and (c+d)
        Returns: (a+b)*v = a*v + b*v for multivector (a+b) and scalar v

        ```python
        >>> from multivectors import x, y, z, w
        >>> (x + y) * (z + w)
        (1.0 * x*z + 1.0 * x*w + 1.0 * y*z + 1.0 * y*w)
        >>> (x + y) * 3
        (3.0 * x + 3.0 * y)
        >>> (x + y) * x
        (1.0 + -1.0 * x*y)

        ```
        """
        if isinstance(other, MultiVector):
            if self.grade is not None and other.grade is not None:
                ((bases1, scalar1),) = self.termdict.items()
                ((bases2, scalar2),) = other.termdict.items()
                (bases, scalar) = condense_bases(
                    (*bases1, *bases2), scalar1 * scalar2)
                return MultiVector({bases: scalar})
            return self.from_terms(*(a * b for b in other.terms
                                     for a in self.terms))
        if isinstance(other, Scalar_t):
            return self.from_terms(*(t * self.from_terms(other) for t in self.terms))
        return NotImplemented

    def __rmul__(self, other: Scalar_a) -> MultiVector:
        """Support multiplying multivectors on the right side of scalars.

        ```python
        >>> from multivectors import x, y
        >>> 3 * (x + y)
        (3.0 * x + 3.0 * y)
        >>> y * (x + y)
        (1.0 + -1.0 * x*y)

        ```
        """
        if isinstance(other, Scalar_t):
            return self.from_terms(*(self.from_terms(other) * t
                                     for t in self.terms))
        return NotImplemented

    def __matmul__(self, other: SOV) -> MultiVector:
        """Get the inner (dot) product of two objects.

        Returns: u@v = (u*v)[abs(u.grade-v.grade)] when grade is defined
        Returns: (a+b)@(c+d)=a@c+a@d+b@c+b@d for multivectors (a+b) and (c+d)
        Returns: (a+b)@v = a@v + b@v for multivector (a+b) and scalar v

        ```python
        >>> from multivectors import x, y, z
        >>> (2*x + 3*y) @ (4*x + 5*y)
        (23.0)
        >>> (2*x*y).dot(3*y*z)
        (0.0)
        >>> (x + y).inner(3)
        (3.0 * x + 3.0 * y)
        >>> (x + y) @ x
        (1.0)

        ```
        """
        if isinstance(other, MultiVector):
            if self.grade is not None and other.grade is not None:
                return (self * other)[abs(self.grade - other.grade)]
            return self.from_terms(*(
                a @ b for b in other.terms for a in self.terms))
        if isinstance(other, Scalar_t):
            return self.from_terms(*(t @ self.from_terms(other)
                                     for t in self.terms))
        return NotImplemented

    dot = inner = __matmul__

    def __rmatmul__(self, other: Scalar_a) -> MultiVector:
        """Support dotting multivectors on the right hand side.

        Returns: v@(a+b) = v@a + v@b for multivector (a+b) and scalar v

        ```python
        >>> from multivectors import x, y
        >>> 3 @ (x + y)
        (3.0 * x + 3.0 * y)
        >>> x @ (x + y)
        (1.0)

        ```
        """
        if isinstance(other, Scalar_t):
            return self.from_terms(*(self.from_terms(other) @ t
                                     for t in self.terms))
        return NotImplemented

    def __truediv__(self, other: SOV) -> MultiVector:
        """Divide two objects.

        Returns: (a+b)/v = a/v + b/v

        ```python
        >>> from multivectors import x, y
        >>> (6*x + 9*y) / 3
        (2.0 * x + 3.0 * y)
        >>> (6*x + 9*y) / (3*x)
        (2.0 + -3.0 * x*y)

        ```
        """
        if not isinstance(other, Scalar_t):
            if not isinstance(other, MultiVector):
                return NotImplemented
            if other.grade is None:
                return NotImplemented
        other = 1 / other
        return MultiVector.from_terms(*(t * other for t in self.terms))

    def __rtruediv__(self, other: Scalar_a) -> MultiVector:
        """Divide a scalar by a multivector. Only defined for blades.

        ```python
        >>> from multivectors import x, y
        >>> 1 / x
        (1.0 * x)
        >>> 2 / (4 * x*y)
        (-0.5 * x*y)

        ```
        """
        if self.grade is None:
            raise TypeError('cannot take inverse of multivector '
                            'with more than one term')
        return other * self / (self * self)._

    def __mod__(self, idxs: Index) -> Scalar_a:
        """Support index swizzling.

        ```python
        >>> from multivectors import x, y
        >>> (x + y) % 0
        1.0
        >>> v = 1 + 2*y + 3*x*y
        >>> v % ()
        1.0
        >>> v % 1
        2.0
        >>> v % (0, 1)
        3.0
        >>> v % 2
        0.0

        ```
        """
        return self.termdict.get(tuple(idxs_to_idxs(idxs)), Scalar_f())

    def __pow__(self, other: int) -> MultiVector:
        """A multivector raised to an integer power.
        V ** n = V * V * V * ... * V n times.
        V ** -n = 1 / (V ** n)

        ```python
        >>> from multivectors import x, y
        >>> (x + y) ** 3
        (2.0 * x + 2.0 * y)
        >>> (2 * x*y) ** -5
        (-0.03125 * x*y)

        """
        if not isinstance(other, int):
            return NotImplemented
        result = MultiVector({(): Scalar_f('1')})
        for _ in range(abs(other)):
            result *= self
        if other < 0:
            return Scalar_f('1') / result
        return result

    def __rpow__(self, other: Scalar_a) -> MultiVector:
        """A real number raised to a multivector power.
        x ** V = Product of x ** V_i = Product of e ** (V_i ln x)
        = Product of (cos(a_i ln x) + sin(a_i ln x) * Ii)

        ```python
        >>> from math import e, pi
        >>> from multivectors import x, y
        >>> round(e ** (pi/4 * (x + y)), 2)
        (0.5 + 0.5 * x + 0.5 * y + 0.5 * x*y)
        >>> round(e ** (pi * (x + y)), 2)
        (1.0)

        ```
        """
        if not isinstance(other, Scalar_t):
            return NotImplemented
        result = MultiVector({(): Scalar_f('1')})
        for t in self.terms:
            ((bases, scalar),) = t.termdict.items()
            theta = scalar * Scalar_f(math.log(other))
            result *= Scalar_f(math.cos(theta)) + MultiVector({
                bases: Scalar_f(math.sin(theta))})
        return result

    def __xor__(self, other: SOV) -> MultiVector:
        """Get the outer (wedge) product of two objects.
        WARNING: Operator precedence puts ^ after +!
        Make sure to put outer products in parentheses, like this:
        ``u * v == u @ v + (u ^ v)``

        Returns: u^v = (u*v)[u.grade+v.grade] when grade is defined
        Returns: (a+b)^(c+d)=(a^c)+(a^d)+(b^c)+(b^d) for (a+b), (c+d)
        Returns: (a+b)^v = (a^v) + (b^v) for multivector (a+b) and scalar v

        ```python
        >>> from multivectors import x, y, z
        >>> (2*x + 3*y) ^ (4*x + 5*y)
        (-2.0 * x*y)
        >>> (2*x*y).wedge(3*y*z)
        (0.0)
        >>> (x + y).outer(3)
        (3.0 * x + 3.0 * y)
        >>> (x + y) ^ x
        (-1.0 * x*y)

        ```
        """
        if isinstance(other, MultiVector):
            if self.grade is not None and other.grade is not None:
                return (self * other)[self.grade + other.grade]
            return self.from_terms(*(
                a ^ b for b in other.terms for a in self.terms))
        if isinstance(other, Scalar_t):
            return self.from_terms(*(t ^ self.from_terms(other)
                                     for t in self.terms))
        return NotImplemented

    wedge = outer = __xor__

    def __rxor__(self, other: Scalar_a) -> MultiVector:
        """Support wedging multivectors on the right hand side.

        Returns: v@(a+b) = v@a + v@b for multivector (a+b) and simple v

        ```python
        >>> from multivectors import x, y
        >>> 3 ^ (x + y)
        (3.0 * x + 3.0 * y)
        >>> x ^ (x + y)
        (1.0 * x*y)

        ```
        """
        if isinstance(other, Scalar_t):
            return self.from_terms(*(self.from_terms(other) ^ t
                                     for t in self.terms))
        return NotImplemented

    # Unary operators

    def __neg__(self) -> MultiVector:
        """The negation of a multivector is the negation of all its terms.

        ```python
        >>> from multivectors import x, y, z
        >>> -(2*x + 3*y*z)
        (-2.0 * x + -3.0 * y*z)

        ```
        """
        return MultiVector({b: -s for b, s in self.termdict.items()})

    def __pos__(self) -> MultiVector:
        """A normalized multivector is one scaled down by its magnitude.

        ```python
        >>> from multivectors import x, y, z, w
        >>> +(x + y + z + w)
        (0.5 * x + 0.5 * y + 0.5 * z + 0.5 * w)
        >>> round((x + y).normalize(), 3)
        (0.707 * x + 0.707 * y)

        ```
        """
        return self / abs(self)

    normalize = __pos__

    def __abs__(self) -> Scalar_a:
        """The magnitude of a multivector is the square root of
        the sum of the squares of its terms.

        ```python
        >>> from multivectors import x, y, z, w
        >>> abs(x + y + z + w)
        2.0
        >>> (2 ** 1.5 * x + 2 ** 1.5 * y).magnitude()
        4.0

        ```
        """
        return Scalar_f(sum((t * t)._ for t in self.terms)) ** Scalar_f('0.5')

    magnitude = __abs__

    def __invert__(self) -> MultiVector:
        """The conjugate of a multivector is the negation of all terms
        besides the real component.

        ```python
        >>> from multivectors import x, y
        >>> ~(1 + 2*x + 3*x*y)
        (1.0 + -2.0 * x + -3.0 * x*y)
        >>> (2 + 3*y).conjugate()
        (2.0 + -3.0 * y)

        ```
        """
        return MultiVector({
            bases: (-scalar if bases else scalar)
            for bases, scalar in self.termdict.items()})

    conjugate = __invert__

    def __complex__(self) -> complex:
        """Convert a scalar to a complex number.

        ```python
        >>> from multivectors import _, x
        >>> complex(2 * _)
        (2+0j)
        >>> complex(3 * x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * x) to complex

        ```
        """
        if self.grade != 0:
            raise TypeError(
                'cannot convert non-scalar blade %r to complex' % self)
        return complex(self._)

    def __int__(self) -> int:
        """Convert a scalar to an integer.

        ```python
        >>> from multivectors import _, x
        >>> int(2 * _)
        2
        >>> int(3 * x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * x) to int

        ```
        """
        if self.grade != 0:
            raise TypeError(
                'cannot convert non-scalar blade %r to int' % self)
        return int(self._)

    def __float__(self) -> float:
        """Convert a scalar to a float.

        ```python
        >>> from multivectors import _, x
        >>> float(2 * _)
        2.0
        >>> float(3 * x)
        Traceback (most recent call last):
            ...
        TypeError: cannot convert non-scalar blade (3.0 * x) to float

        ```
        """
        if self.grade != 0:
            raise TypeError(
                'cannot convert non-scalar blade %r to float' % self)
        return float(self._)

    def __round__(self, ndigits: int = None) -> MultiVector:
        """Round the scalars of each component term of a multivector.

        ```python
        >>> from multivectors import x, y
        >>> round(1.7 * x + 1.2 * y)
        (2.0 * x + 1.0 * y)
        >>> round(0.15 * x + 0.05 * y, 1)
        (0.1 * x + 0.1 * y)

        ```
        """
        return MultiVector({b: round(s, ndigits)
                            for b, s in self.termdict.items()})

    def __trunc__(self) -> MultiVector:
        """Truncate the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> from multivectors import x, y
        >>> math.trunc(1.7 * x + 1.2 * y)
        (1.0 * x + 1.0 * y)
        >>> math.trunc(-1.7 * x - 1.2 * y)
        (-1.0 * x + -1.0 * y)

        ```
        """
        return MultiVector({b: math.trunc(s)
                            for b, s in self.termdict.items()})

    def __floor__(self) -> MultiVector:
        """Floor the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> from multivectors import x, y
        >>> math.floor(1.7 * x + 1.2 * y)
        (1.0 * x + 1.0 * y)
        >>> math.floor(-1.7 * x - 1.2 * y)
        (-2.0 * x + -2.0 * y)

        ```
        """
        return MultiVector({b: math.floor(s)
                            for b, s in self.termdict.items()})

    def __ceil__(self) -> MultiVector:
        """Ceiling the scalars of each component term of a multivector.

        ```python
        >>> import math
        >>> from multivectors import x, y
        >>> math.ceil(1.7 * x + 1.2 * y)
        (2.0 * x + 2.0 * y)
        >>> math.ceil(-1.7 * x - 1.2 * y)
        (-1.0 * x + -1.0 * y)

        ```
        """
        return MultiVector({b: math.ceil(s)
                            for b, s in self.termdict.items()})

    # Actual methods

    def rotate(self, angle: Scalar_a, plane: MultiVector) -> MultiVector:
        """Rotate this multivector by angle in rads around the blade plane.

        angle: Angle to rotate by, in radians.
        plane: Blade representing basis plane to rotate through.

        ```python
        >>> from math import radians
        >>> from multivectors import x, y, z, w
        >>> round((3*x + 2*y + 4*z).rotate(
        ...     radians(90), x*y), 2)
        (-2.0 * x + 3.0 * y + 4.0 * z)
        >>> round((3*x + 2*y + 4*z + 5*w).rotate(
        ...     radians(90), x*y*z), 2)
        (3.0 * x + 2.0 * y + 4.0 * z + -5.0 * x*y*z*w)

        ```
        """
        if plane.grade is None or abs(plane * plane) != 1:
            raise TypeError('%s is not a basis plane' % plane)
        R = math.e ** (plane * angle / 2)
        return ~R * self * R

    def angle_to(self, other: MultiVector) -> Scalar_a:
        """Get the angle between this multivector and another.

        ```python
        >>> from math import degrees
        >>> from multivectors import x, y, z, w
        >>> math.degrees((x + y).angle_to(x - y))
        90.0
        >>> round(math.degrees((x + y + z + w).angle_to(x - y - z - w)), 2)
        120.0

        ```
        """
        return math.acos((self @ other) / (abs(self) * abs(other)))

Index = Union[int, Iterable[int], slice]
Scalar_a = Real # annotation
Scalar_t = Real # isinstance()d type
Scalar_f = float # factory
TermDict = Dict[Tuple[int, ...], Scalar_a]
SOV = Union[Scalar_a, MultiVector]

_blades: Dict[str, MultiVector] = {}

def __getattr__(name: str) -> MultiVector:
    """Support module level swizzling of basis vectors.
    Unlike Blade.__getattr__, this rejects invalid characters in the name.
    """
    return _blades.setdefault(name, MultiVector(
        {tuple(names_to_idxs(name, True)): Scalar_f('1')}))

def set_scalar_factory(
    instancecheck: Union[Type, Tuple[Type, ...]],
    factory: Type,
    annotation: Type = None,
):
    """Set the factory for scalars.

    This could be e.g. Decimal for fixed-precision scalars,
    or Fraction for rational scalars. The default factory is `float`.

    annotation: Function and variable annotations use this.
    instancecheck: Used in isinstance(), may be type or tuple of types.
        Should include `numbers.Real` if you want to preserve arithmetic
        with literals like `2 * x`.
    factory: See below.

    The `factory` provided must be a constructor that implements the following:
    - `factory()` = the 0 value of the class
    - `factory(f: float)` = the float value converted to the class
    - `factory(s: str)` = the string value parsed to the class

    This will not update any existing multivectors (including module-getattred
    ones) so make sure to re-get them after calling this function!
    """
    global Scalar_a, Scalar_f, Scalar_t
    Scalar_a = annotation or Union[instancecheck]
    Scalar_f = factory
    Scalar_t = instancecheck
    _blades.clear()
