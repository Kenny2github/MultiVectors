"""
Compute multivectors in arbitrary-dimensional space.

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
"""
from __future__ import annotations
import math
from numbers import Real
from typing import Dict, Iterable, List, Optional, Tuple, Union

__all__ = [
    'MultiVector',
    '_',
    'x',
    'y',
    'z',
    'w'
]

__version__ = '1.2.1'

NAMES = 'xyzw'

def merge(arr: List[int], left: List[int], right: List[int]) -> int:
    """Perform a merge of sorted lists and count the swaps.

    Args:
        arr: The list to write into.
        left: The left sorted list.
        right: The right sorted list.

    Returns:
        The number of swaps made when merging.
    """
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

    Args:
        arr: The list to sort.
        copy: If :obj:`True`, (the default) don't modify the original list.

    Returns:
        The number of swaps made when sorting the list.

    Examples:

        >>> count_swaps([1, 3, 2, 5, 4])
        2
        >>> count_swaps([3, 2, 1])
        3
    """
    if copy:
        arr = list(arr)
    if len(arr) < 2:
        return 0
    mid = (len(arr) + 1) // 2
    left = arr[:mid]
    right = arr[mid:]
    return (count_swaps(left, False)
            + count_swaps(right, False)
            + merge(arr, left, right))

def names_to_idxs(name: str,
                  raise_on_invalid_chars: bool = False) -> List[int]:
    """Convert swizzled basis vector names into generalized basis indices.

    Args:
        name: The names to convert.
        raise_on_invalid_chars: If :obj:`True` (default :obj:`False`),
            raise :exc:`AttributeError` if any characters appear in the
            name that are invalid in a basis name.

    Returns:
        A list of basis indexes that the name represents.

    Raises:
        AttributeError: If characters invalid for a basis name appear,
            and ``raise_on_invalid_chars`` is :obj:`True`.

    Examples:

        >>> names_to_idxs('xyw')
        [0, 1, 3]
        >>> names_to_idxs('e_1e2_z')
        [0, 1, 2]
        >>> names_to_idxs('_')
        []
    """
    if name.startswith('__'):
        # fail on magic attributes immediately
        raise AttributeError
    idxs: List[int] = []
    index = 0
    name_len = len(name)
    while index < name_len:
        char = name[index]
        if char in NAMES:
            idxs.append(NAMES.index(char))
            index += 1
        elif char == 'e':
            index += 1
            end = index
            subname = []
            while (
                end < name_len
                and name[end] not in (NAMES + 'e')
                and not subname
            ):
                if '0' <= name[end] <= '9':
                    subname.append(name[end])
                elif raise_on_invalid_chars:
                    raise AttributeError
                end += 1
            try:
                idxs.append(int(''.join(subname)) - 1)
            except ValueError:
                raise AttributeError(f'empty eN notation: {name!r}') from None
            index = end
        elif char == '_':
            index += 1
        elif raise_on_invalid_chars:
            raise AttributeError
        else:
            index += 1
    return idxs

def idxs_to_idxs(idxs: Index) -> List[int]:
    """Convert multiple possible ways to specify multiple indices.

    This is intended to be given the argument to :class:`MultiVector` indexing:
    ``V[0]`` would call ``idxs_to_idxs(0)``, ``V[0, 1]`` would call
    ``idxs_to_idxs((0, 1))``, and ``V[0:1]`` would call
    ``idxs_to_idxs(slice(0, 1, None))``

    Args:
        idxs: The indexes to convert. This can be an integer, for just that
            index; an iterable of integers, for those indexes directly; or a
            slice, for the indexes the slice represents.

    Returns:
        A list of basis indexes, converted from the argument.

    Examples:

        >>> idxs_to_idxs(slice(None, 5, None))
        [0, 1, 2, 3, 4]
        >>> idxs_to_idxs((1, 3, 4))
        [1, 3, 4]
        >>> idxs_to_idxs(1)
        [1]
    """
    if isinstance(idxs, int):
        return [idxs]
    if isinstance(idxs, slice):
        if idxs.stop is None:
            raise TypeError('cannot have infinite bases')
        return list(range(*idxs.indices(idxs.stop)))
    return list(idxs)

def idxs_to_names(idxs: Index, sep: str = '') -> str:
    """Convert indices to a swizzled name combination.

    Args:
        idxs: The basis index(es), as accepted by :func:`idxs_to_idxs`.
        sep: The separator to use between the basis vector names.

    Returns:
        The swizzled names.

    Examples:

        >>> idxs_to_names(slice(None, 5, None))
        'e1e2e3e4e5'
        >>> idxs_to_names((0, 1, 3))
        'xyw'
        >>> idxs_to_names(2)
        'z'
        >>> idxs_to_names((0, 2, 3), sep='*')
        'x*z*w'
    """
    idxs = idxs_to_idxs(idxs)
    if max(idxs) > 3:
        return sep.join(f'e{idx + 1}' for idx in idxs)
    return sep.join(NAMES[idx] for idx in idxs)

def condense_bases(bases: Tuple[int, ...], scalar: float = 1.0) \
        -> Tuple[Tuple[int, ...], float]:
    """Normalize a sequence of bases, modifying the scalar as necessary.

    Args:
        bases: The tuple of basis indices.
        scalar: Real number that will scale the resulting bases.

    Returns:
        A 2-tuple of normalized bases and the modified scalar.

    Examples:

        >>> condense_bases((1, 1, 2, 1, 2), 2.0)
        ((1,), -2.0)
        >>> condense_bases((1, 2, 1, 2), 1.5)
        ((), -1.5)
        >>> condense_bases((2, 1, 3, 2, 3, 3), 1.0)
        ((1, 3), 1.0)
    """
    _bases = list(bases)
    if count_swaps(_bases) % 2 != 0:
        scalar *= -1.0
    _bases.sort()
    _bases = [basis for basis in set(_bases)
              if _bases.count(basis) % 2 != 0]
    return tuple(_bases), scalar

class MultiVector:
    """A linear combination of geometric products of basis vectors.

    The bare constructor is not meant for regular use.
    Use module swizzling or the factory
    :meth:`~MultiVector.from_terms()` instead.

    Basis vector names can be swizzled on instances:

        >>> from multivectors import x, y, z
        >>> (x + y).x
        1.0
        >>> (x*y + z).xy
        1.0
        >>> (x + y).e3
        0.0

    And indices can be combined:

        >>> (x + y) % 0
        1.0
        >>> (x*y + z) % (0, 1)
        1.0
        >>> (x + y) % 2
        0.0
    """

    termdict: TermDict

    @property
    def grade(self) -> Union[int, None]:
        """The grade of this blade.

        Returns:
            The number of different bases this blade consists of,
            or :obj:`None` if this multivector is not a blade (one term).

        Examples:

            >>> from multivectors import x, y, z
            >>> (x + y).grade # not a blade
            >>> (z * 2 + z).grade
            1
            >>> (x*y*z).grade
            3
        """
        if len(self.termdict) != 1:
            return None
        (bases,) = self.termdict.keys()
        return len(bases)

    @property
    def terms(self) -> Tuple[MultiVector, ...]:
        """Get a sequence of blades comprising this multivector.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> (x + y).terms
            ((1.0 * x), (1.0 * y))
            >>> ((x + y) * (z + w)).terms
            ((1.0 * x*z), (1.0 * x*w), (1.0 * y*z), (1.0 * y*w))
        """
        items = sorted(self.termdict.items(),
                       key=lambda item: (len(item[0]), item[0]))
        return tuple(MultiVector({key: value}) for key, value in items)

    def __init__(self, termdict: TermDict):
        """:meta private:"""
        self.termdict = {bases: float(scalar)
                         for bases, scalar in termdict.items()
                         if scalar != 0.0} or {(): 0.0}

    @classmethod
    def from_terms(cls, *terms: SOV) -> MultiVector:
        """Create a multivector by summing a sequence of terms.

        Args:
            *terms: The terms. If you have an iterable of terms, use
                (e.g.) ``from_terms(*terms)``

        Returns:
            A multivector.

        Examples:

            >>> from multivectors import x, y, z
            >>> MultiVector.from_terms(x, y)
            (1.0 * x + 1.0 * y)
            >>> MultiVector.from_terms()
            (0.0)
            >>> MultiVector.from_terms(z)
            (1.0 * z)
            >>> MultiVector.from_terms(2 * x, x)
            (3.0 * x)
        """
        termseq = [
            t for term in terms for t in
            (term if isinstance(term, MultiVector)
             else MultiVector({(): term})).termdict.items()
        ]
        termsdict: Dict[Tuple[int, ...], List[float]] = {}
        for bases, scalar in termseq:
            termsdict.setdefault(bases, []).append(scalar)
        termdict = {bases: math.fsum(scalars)
                    for bases, scalars in termsdict.items()}
        # discard zeros
        termdict = {bases: scalar for bases, scalar in termdict.items()
                    if scalar != 0.0}
        return cls(termdict)

    @classmethod
    def scalar(cls, num) -> MultiVector:
        """Create a MultiVector representing a scalar.

        Args:
            num: Any object that can be :func:`float()`-d.

        Returns:
            A multivector with only a scalar part of ``num``.

        Examples:

            >>> from multivectors import MultiVector
            >>> MultiVector.scalar('0')
            (0.0)
            >>> MultiVector.scalar('1.2')
            (1.2)
        """
        return cls({(): float(num)})

    def __getattr__(self, name: str) -> float:
        """Support basis name swizzling."""
        return self.termdict.get(tuple(names_to_idxs(name)), 0.0)

    def __getitem__(self, grades: Index) -> MultiVector:
        """The choose operator - returns the sum of all blades of grade k.

        Args:
            grades: The grade(s) to choose.

        Examples:

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
            >>> (1 + 2*x + 3*x*y)[1:]
            (2.0 * x + 3.0 * x*y)
        """
        if isinstance(grades, slice):
            stop = (max(map(len, self.termdict.keys())) + 1) \
                if grades.stop is None else grades.stop
            grades = slice(grades.start, stop, grades.step)
        grades = set(idxs_to_idxs(grades))
        return MultiVector({
            bases: scalar for bases, scalar in self.termdict.items()
            if len(bases) in grades})

    choose = __getitem__

    def __repr__(self) -> str:
        """Return a representation of this multivector.
        Depending on the global namespace, this may be :func:`eval()`-able.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> repr(x)
            '(1.0 * x)'
            >>> repr(x + y)
            '(1.0 * x + 1.0 * y)'
            >>> repr(y*z - x*w)
            '(-1.0 * x*w + 1.0 * y*z)'
        """
        if self.grade is not None:
            ((bases, scalar),) = self.termdict.items()
            if bases == ():
                return f'({scalar!r})'
            names = idxs_to_names(bases, '*')
            return f'({scalar!r} * {names!s})'
        return '(' + ' + '.join(repr(term).strip('()')
                                for term in self.terms) + ')'

    def __str__(self) -> str:
        """Return a representation of this multivector suited for showing.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> str(x)
            '1.00x'
            >>> str(x + y)
            '(1.00x + 1.00y)'
            >>> str(y*z - x*w)
            '(-1.00xw + 1.00yz)'
            >>> print(1 + x + x*y)
            (1.00 + 1.00x + 1.00xy)
        """
        return f'{self:.2f}'

    def __format__(self, spec: str) -> str:
        """Return a representation of this multivector suited for formatting.

        Args:
            spec: The format spec (forwarded to the underlying :class:`float`s)

        Examples:

            >>> from multivectors import x, y, z, w
            >>> f'{x:.3f}'
            '1.000x'
            >>> V = x + 2 * y + z
            >>> '{:.3f}'.format(V)
            '(1.000x + 2.000y + 1.000z)'
            >>> V = y*z - x*w
            >>> f'{V:.1f}'
            '(-1.0xw + 1.0yz)'
        """
        if self.grade is not None:
            ((bases, scalar),) = self.termdict.items()
            if bases == ():
                return format(scalar, spec)
            names = idxs_to_names(bases)
            return format(scalar, spec) + str(names)
        return '(' + ' + '.join(format(term, spec) for term in self.terms) + ')'

    # Relational operators

    def __eq__(self, other: SOV) -> bool:
        """Compare equality of two objects.

        Returns:
            :obj:`True` if all terms of this multivector are equal to the
            ``other``; :obj:`True` if this multivector is scalar and equals
            the ``other``; or :obj:`False` for all other cases or types.

        Examples:

            >>> from multivectors import x, y
            >>> x + y == y + x
            True
            >>> x + 2*y == 2*x + y
            False
        """
        if not isinstance(other, MultiVector):
            if self.grade != 0:
                return False # this MV is not a scalar
            return self._ == other
        return self.termdict == other.termdict

    def __ne__(self, other: SOV) -> bool:
        """Compare inequality of two objects.

        Returns:
            :obj:`False` if all terms of this multivector are equal to the
            ``other``; :obj:`False` if this multivector is scalar and equals
            the ``other``; or :obj:`True` for all other cases or types.

        Examples:

            >>> from multivectors import x, y
            >>> x + y != y + x
            False
            >>> x + 2*y != 2*x + y
            True
        """
        return not (self == other)

    def __lt__(self, other: float) -> bool:
        """Compare this blade less than an object.

        Returns:
            :obj:`True` if this is a scalar blade less than the scalar;
            :obj:`False` if this is a scalar blade not less than the scalar;
            or :obj:`NotImplemented` for all other types.

        Examples:

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
        """
        if not isinstance(other, Real) or self.grade != 0:
            return NotImplemented
        return self._ < other

    def __gt__(self, other: float) -> bool:
        """Compare this blade greater than an object.

        Returns:
            :obj:`True` if this is a scalar blade greater than the scalar;
            :obj:`False` if this is a scalar blade not greater than the scalar;
            or :obj:`NotImplemented` for all other types.

        Examples:

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
        """
        if not isinstance(other, Real) or self.grade != 0:
            return NotImplemented
        return self._ > other

    def __le__(self, other: float) -> bool:
        """Compare this blade less than or equal to an object.

        Returns:
            :obj:`True` if this is a scalar blade less than or equal to the
            scalar; :obj:`False` if this is a scalar blade greater than the
            scalar; or :obj:`NotImplemented` for all other types.

        Examples:

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
        """
        if not isinstance(other, Real) or self.grade != 0:
            return NotImplemented
        return self._ <= other

    def __ge__(self, other: float) -> bool:
        """Compare this blade greater than or equal to an object.

        Returns:
            :obj:`True` if this is a scalar blade greater than or equal to the
            scalar; :obj:`False` if this is a scalar blade less than the
            scalar; or :obj:`NotImplemented` for all other types.

        Examples:

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
        """
        if not isinstance(other, Real) or self.grade != 0:
            return NotImplemented
        return self._ >= other

    # Binary operators

    def __add__(self, other: SOV) -> MultiVector:
        """Add a multivector and another object.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> (x + z) + (y + w)
            (1.0 * x + 1.0 * y + 1.0 * z + 1.0 * w)
            >>> (x + z) + y
            (1.0 * x + 1.0 * y + 1.0 * z)
            >>> (x + y) + 1
            (1.0 + 1.0 * x + 1.0 * y)
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*self.terms, *other.terms)
        if isinstance(other, Real):
            return self.from_terms(*self.terms, other)
        return NotImplemented

    def __radd__(self, other: float) -> MultiVector:
        """Support adding multivectors on the right side of objects.

        Examples:

            >>> from multivectors import x, y, z
            >>> 1 + (x + z)
            (1.0 + 1.0 * x + 1.0 * z)
            >>> x + (y + z)
            (1.0 * x + 1.0 * y + 1.0 * z)
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*other.terms, *self.terms)
        if isinstance(other, Real):
            return self.from_terms(other, *self.terms)
        return NotImplemented

    def __sub__(self, other: SOV) -> MultiVector:
        """Subtracting is adding the negation.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> (x + z) - (y + w)
            (1.0 * x + -1.0 * y + 1.0 * z + -1.0 * w)
            >>> (x + z) - y
            (1.0 * x + -1.0 * y + 1.0 * z)
        """
        return self + (-other)

    def __rsub__(self, other: float) -> MultiVector:
        """Support subtracting multivectors from objects.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> 1 - (x + z)
            (1.0 + -1.0 * x + -1.0 * z)
            >>> x - (y + z)
            (1.0 * x + -1.0 * y + -1.0 * z)
        """
        return other + (-self)

    def __mul__(self, other: SOV) -> MultiVector:
        """Multiply a multivector and another object.

        Returns:
            ``(a + b) * (c + d) = a*c + a*d + b*c + b*d``
            for multivectors ``(a + b)`` and ``(c + d)``.

        Returns:
            ``(a + b) * v = a*v + b*v``
            for multivector ``(a + b)`` and scalar ``v``.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> (x + y) * (z + w)
            (1.0 * x*z + 1.0 * x*w + 1.0 * y*z + 1.0 * y*w)
            >>> (x + y) * 3
            (3.0 * x + 3.0 * y)
            >>> (x + y) * x
            (1.0 + -1.0 * x*y)
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
        if isinstance(other, Real):
            return self.from_terms(*(term * self.from_terms(other)
                                     for term in self.terms))
        return NotImplemented

    def __rmul__(self, other: float) -> MultiVector:
        """Support multiplying multivectors on the right side of scalars.

        Examples:

            >>> from multivectors import x, y
            >>> 3 * (x + y)
            (3.0 * x + 3.0 * y)
            >>> y * (x + y)
            (1.0 + -1.0 * x*y)
        """
        if isinstance(other, Real):
            return self.from_terms(*(self.from_terms(other) * term
                                     for term in self.terms))
        return NotImplemented

    def __matmul__(self, other: SOV) -> MultiVector:
        """Get the inner (dot) product of two objects.

        Returns:
            ``u @ v = (u * v)[abs(u.grade - v.grade)]``
            when :attr:`~MultiVector.grade` is defined.

        Returns:
            ``(a + b) @ (c + d) = a@c + a@d + b@c + b@d``
            for multivectors ``(a + b)`` and ``(c + d)``

        Returns:
            ``(a + b) @ v = a@v + b@v``
            for multivector ``(a + b)`` and scalar ``v``

        Examples:

            >>> from multivectors import x, y, z
            >>> (2*x + 3*y) @ (4*x + 5*y)
            (23.0)
            >>> (2*x*y).dot(3*y*z)
            (0.0)
            >>> (x + y).inner(3)
            (3.0 * x + 3.0 * y)
            >>> (x + y) @ x
            (1.0)
        """
        if isinstance(other, MultiVector):
            if self.grade is not None and other.grade is not None:
                return (self * other)[abs(self.grade - other.grade)]
            return self.from_terms(*(
                a @ b for b in other.terms for a in self.terms))
        if isinstance(other, Real):
            return self.from_terms(*(term @ self.from_terms(other)
                                     for term in self.terms))
        return NotImplemented

    dot = inner = __matmul__

    def __rmatmul__(self, other: float) -> MultiVector:
        """Support dotting multivectors on the right hand side.

        Returns:
            ``v @ (a + b) = v@a + v@b``
            for multivector ``(a + b)`` and scalar ``v``

        Examples:

            >>> from multivectors import x, y
            >>> 3 @ (x + y)
            (3.0 * x + 3.0 * y)
            >>> x @ (x + y)
            (1.0)
        """
        if isinstance(other, Real):
            return self.from_terms(*(self.from_terms(other) @ term
                                     for term in self.terms))
        return NotImplemented

    def __truediv__(self, other: SOV) -> MultiVector:
        """Divide two objects.

        Returns:
            ``(a + b) / v = a/v + b/v``

        Examples:

            >>> from multivectors import x, y
            >>> (6*x + 9*y) / 3
            (2.0 * x + 3.0 * y)
            >>> (6*x + 9*y) / (3*x)
            (2.0 + -3.0 * x*y)
        """
        if not isinstance(other, Real):
            if not isinstance(other, MultiVector):
                return NotImplemented
            if other.grade is None:
                return NotImplemented
        other = 1 / other
        return MultiVector.from_terms(*(term * other for term in self.terms))

    def __rtruediv__(self, other: float) -> MultiVector:
        """Divide a scalar by a multivector. Only defined for blades.

        Examples:

            >>> from multivectors import x, y
            >>> 1 / x
            (1.0 * x)
            >>> 2 / (4 * x*y)
            (-0.5 * x*y)
        """
        if self.grade is None:
            raise TypeError('cannot take inverse of multivector '
                            'with more than one term')
        return other * self / (self * self)._

    def __mod__(self, idxs: Index) -> float:
        """Support index swizzling.

        Examples:

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
        """
        return self.termdict.get(tuple(idxs_to_idxs(idxs)), 0.0)

    def __pow__(self, other: int) -> MultiVector:
        """A multivector raised to an integer power.

        ``V ** n = V * V * V * ... * V``, ``n`` times.
        ``V ** -n = 1 / (V ** n)``

        Examples:

            >>> from multivectors import x, y
            >>> (x + y) ** 3
            (2.0 * x + 2.0 * y)
            >>> (2 * x*y) ** -5
            (-0.03125 * x*y)
        """
        if not isinstance(other, int):
            return NotImplemented
        result = self.scalar(1)
        for _ in range(abs(other)):
            result *= self
        if other < 0:
            return 1.0 / result
        return result

    def __rpow__(self, other: float) -> MultiVector:
        """A real number raised to a multivector power.

        ``x ** V = e ** ln(x ** V) = e ** (V ln x)``

        Examples:

            >>> from multivectors import x, y
            >>> round(2 ** (x + y), 2)
            (1.52 + 0.81 * x + 0.81 * y)
        """
        if not isinstance(other, Real):
            return NotImplemented
        return (self * math.log(other)).exp()

    def exp(self) -> MultiVector:
        r"""e raised to this multivector power.

        .. math::

            e^V = \exp(V) = \sum_{n=0}^{\infty} \frac{V^n}{n!}

        Examples:

            >>> from math import pi, sqrt
            >>> from multivectors import x, y
            >>> # 45-degree rotation through xy-plane
            >>> # results in (1+xy)/sqrt(2)
            >>> (pi/4 * x*y).exp() * sqrt(2)
            (1.0 + 1.0 * x*y)
            >>> round((pi * x*y).exp(), 14)
            (-1.0)
        """
        # since we have to sequentially sum terms anyway, we might
        # as well repeatedly multiply for the integer exponent
        # and sequentially multiply for the factorial
        i = current_factorial = 1
        result = self.scalar(1)
        current_exponent = self.scalar(1)
        last_result = self.scalar(0)
        while last_result != result:
            last_result = result
            current_exponent *= self
            current_factorial *= i
            i += 1
            result += current_exponent / current_factorial
        return last_result

    def __xor__(self, other: SOV) -> MultiVector:
        """Get the outer (wedge) product of two objects.

        .. warning::

            Operator precedence puts ``^`` after ``+``!
            Make sure to put outer products in parentheses, like this:
            ``u * v == u @ v + (u ^ v)``

        Returns:
            ``u ^ v = (u * v)[u.grade + v.grade]``
            when :attr:`~MultiVector.grade` is defined

        Returns:
            ``(a + b) ^ (c + d) = (a^c) + (a^d) + (b^c) + (b^d)``
            for multivector ``(a + b)`` and ``(c + d)``

        Returns:
            ``(a + b) ^ v = (a^v) + (b^v)``
            for multivector ``(a + b)`` and scalar ``v``

        Examples:

            >>> from multivectors import x, y, z
            >>> (2*x + 3*y) ^ (4*x + 5*y)
            (-2.0 * x*y)
            >>> (2*x*y).wedge(3*y*z)
            (0.0)
            >>> (x + y).outer(3)
            (3.0 * x + 3.0 * y)
            >>> (x + y) ^ x
            (-1.0 * x*y)
        """
        if isinstance(other, MultiVector):
            if self.grade is not None and other.grade is not None:
                return (self * other)[self.grade + other.grade]
            return self.from_terms(*(
                a ^ b for b in other.terms for a in self.terms))
        if isinstance(other, Real):
            return self.from_terms(*(term ^ self.from_terms(other)
                                     for term in self.terms))
        return NotImplemented

    wedge = outer = __xor__

    def __rxor__(self, other: float) -> MultiVector:
        """Support wedging multivectors on the right hand side.

        Returns:
            ``v ^ (a + b) = (v^a) + (v^b)``
            for multivector ``(a + b)`` and simple ``v``

        Examples:

            >>> from multivectors import x, y
            >>> 3 ^ (x + y)
            (3.0 * x + 3.0 * y)
            >>> x ^ (x + y)
            (1.0 * x*y)
        """
        if isinstance(other, Real):
            return self.from_terms(*(self.from_terms(other) ^ term
                                     for term in self.terms))
        return NotImplemented

    # Unary operators

    def __neg__(self) -> MultiVector:
        """The negation of a multivector is the negation of all its terms.

        Examples:

            >>> from multivectors import x, y, z
            >>> -(2*x + 3*y*z)
            (-2.0 * x + -3.0 * y*z)
        """
        return MultiVector({bases: -scalar for bases, scalar
                            in self.termdict.items()})

    def __pos__(self) -> MultiVector:
        """A normalized multivector is one scaled down by its magnitude.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> +(x + y + z + w)
            (0.5 * x + 0.5 * y + 0.5 * z + 0.5 * w)
            >>> round((x + y).normalize(), 3)
            (0.707 * x + 0.707 * y)
        """
        return self / abs(self)

    normalize = __pos__

    def __abs__(self) -> float:
        """The magnitude of a multivector is the square root of
        the sum of the squares of its terms.

        Examples:

            >>> from multivectors import x, y, z, w
            >>> abs(x + y + z + w)
            2.0
            >>> (2 ** 1.5 * x + 2 ** 1.5 * y).magnitude()
            4.0
        """
        return math.fsum((term * term)._ for term in self.terms) ** 0.5

    magnitude = __abs__

    def __invert__(self) -> MultiVector:
        """The conjugate of a multivector is the negation of all terms
        besides the real component.

        Examples:

            >>> from multivectors import x, y
            >>> ~(1 + 2*x + 3*x*y)
            (1.0 + -2.0 * x + -3.0 * x*y)
            >>> (2 + 3*y).conjugate()
            (2.0 + -3.0 * y)
        """
        return MultiVector({
            bases: (-scalar if bases else scalar)
            for bases, scalar in self.termdict.items()})

    conjugate = __invert__

    def __complex__(self) -> complex:
        """Convert a scalar to a complex number.

        Examples:

            >>> from multivectors import _, x
            >>> complex(2 * _)
            (2+0j)
            >>> complex(3 * x)
            Traceback (most recent call last):
                ...
            TypeError: cannot convert non-scalar blade (3.0 * x) to complex
        """
        if self.grade != 0:
            raise TypeError(
                f'cannot convert non-scalar blade {self!r} to complex')
        return complex(self._)

    def __int__(self) -> int:
        """Convert a scalar to an integer.

        Examples:

            >>> from multivectors import _, x
            >>> int(2 * _)
            2
            >>> int(3 * x)
            Traceback (most recent call last):
                ...
            TypeError: cannot convert non-scalar blade (3.0 * x) to int
        """
        if self.grade != 0:
            raise TypeError(
                f'cannot convert non-scalar blade {self!r} to int')
        return int(self._)

    def __float__(self) -> float:
        """Convert a scalar to a float.

        Examples:

            >>> from multivectors import _, x
            >>> float(2 * _)
            2.0
            >>> float(3 * x)
            Traceback (most recent call last):
                ...
            TypeError: cannot convert non-scalar blade (3.0 * x) to float
        """
        if self.grade != 0:
            raise TypeError(
                f'cannot convert non-scalar blade {self!r} to float')
        return float(self._)

    def __round__(self, ndigits: Optional[int] = None) -> MultiVector:
        """Round the scalars of each component term of a multivector.

        Examples:

            >>> from multivectors import x, y
            >>> round(1.7 * x + 1.2 * y)
            (2.0 * x + 1.0 * y)
            >>> round(0.15 * x + 0.05 * y, 1)
            (0.1 * x + 0.1 * y)
        """
        return MultiVector({bases: round(scalar, ndigits)
                            for bases, scalar in self.termdict.items()})

    def __trunc__(self) -> MultiVector:
        """Truncate the scalars of each component term of a multivector.

        Examples:

            >>> import math
            >>> from multivectors import x, y
            >>> math.trunc(1.7 * x + 1.2 * y)
            (1.0 * x + 1.0 * y)
            >>> math.trunc(-1.7 * x - 1.2 * y)
            (-1.0 * x + -1.0 * y)
        """
        return MultiVector({bases: math.trunc(scalar)
                            for bases, scalar in self.termdict.items()})

    def __floor__(self) -> MultiVector:
        """Floor the scalars of each component term of a multivector.

        Examples:

            >>> import math
            >>> from multivectors import x, y
            >>> math.floor(1.7 * x + 1.2 * y)
            (1.0 * x + 1.0 * y)
            >>> math.floor(-1.7 * x - 1.2 * y)
            (-2.0 * x + -2.0 * y)
        """
        return MultiVector({bases: math.floor(scalar)
                            for bases, scalar in self.termdict.items()})

    def __ceil__(self) -> MultiVector:
        """Ceiling the scalars of each component term of a multivector.

        Examples:

            >>> import math
            >>> from multivectors import x, y
            >>> math.ceil(1.7 * x + 1.2 * y)
            (2.0 * x + 2.0 * y)
            >>> math.ceil(-1.7 * x - 1.2 * y)
            (-1.0 * x + -1.0 * y)
        """
        return MultiVector({bases: math.ceil(scalar)
                            for bases, scalar in self.termdict.items()})

    # Actual methods

    def rotate(self, angle: float, plane: MultiVector) -> MultiVector:
        """Rotate this multivector by angle in rads around the blade plane.

        Args:
            angle: Angle to rotate by, in radians.
            plane: Blade representing basis plane to rotate through.

        Returns:
            Rotated multivector.

        Examples:

            >>> from math import radians
            >>> from multivectors import x, y, z, w
            >>> round((3*x + 2*y + 4*z).rotate(
            ...     radians(90), x*y), 2)
            (-2.0 * x + 3.0 * y + 4.0 * z)
            >>> round((3*x + 2*y + 4*z + 5*w).rotate(
            ...     radians(90), x*y*z), 2)
            (3.0 * x + 2.0 * y + 4.0 * z + -5.0 * x*y*z*w)
        """
        if plane.grade is None or abs(plane * plane) != 1:
            raise TypeError(f'{plane!s} is not a basis plane')
        R = (plane * angle / 2).exp()
        return ~R * self * R

    def angle_to(self, other: MultiVector) -> float:
        """Get the angle between this multivector and another.

        Examples:

            >>> from math import degrees
            >>> from multivectors import x, y, z, w
            >>> math.degrees((x + y).angle_to(x - y))
            90.0
            >>> round(math.degrees((x + y + z + w).angle_to(x - y - z - w)), 2)
            120.0
        """
        return math.acos((self @ other) / (abs(self) * abs(other)))

# type aliases

Index = Union[int, Iterable[int], slice]
TermDict = Dict[Tuple[int, ...], float]
SOV = Union[float, MultiVector]

# module level swizzling

_blades: Dict[str, MultiVector] = {}

def __getattr__(name: str) -> MultiVector:
    """Support module level swizzling of basis vectors.

    Unlike :meth:`MultiVector.__getattr__`,
    this rejects invalid characters in the name.
    """
    return _blades.setdefault(name, MultiVector(
        {tuple(names_to_idxs(name, True)): 1.0}))

_: MultiVector
x: MultiVector
y: MultiVector
z: MultiVector
w: MultiVector
