from __future__ import annotations
import itertools
from typing import Dict, Iterable, List, MutableMapping, Tuple, Union
from itertools import combinations

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
    """Count the number of swaps needed to sort a list."""
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

class Simple:
    """An object with a single grade."""
    bases: List[int]

    @property
    def grade(self) -> int:
        """The grade of this blade."""
        return len(self.bases)

class Scalar(float, Simple):
    """A simple with no bases -- a scalar."""

    def __new__(cls, x):
        if x is NotImplemented:
            return NotImplemented
        return super().__new__(cls, x)

    @property
    def bases(self) -> List[int]:
        return []

    def __add__(self, other: float) -> Scalar:
        return Scalar(super().__add__(other))

    def __sub__(self, other: float) -> Scalar:
        return Scalar(super().__sub__(other))

    def __mul__(self, other: float) -> Scalar:
        return Scalar(super().__mul__(other))

    def __truediv__(self, other: float) -> Scalar:
        return Scalar(super().__truediv__(other))

    def __floordiv__(self, other: float) -> Scalar:
        return Scalar(super().__floordiv__(other))

    def __mod__(self, other: float) -> Scalar:
        return Scalar(super().__mod__(other))

    def __divmod__(self, other: float) -> Scalar:
        return Scalar(super().__divmod__(other))

    def __pow__(self, other: float) -> Scalar:
        return Scalar(super().__pow__(other))

class VectorN:
    """An n-dimensional single-vector."""

    components: List[Scalar]
    names = 'xyzw'

    @property
    def dimension(self) -> int:
        """The number of components in this vector."""
        return len(self.components)

    def __init__(self, *components: Scalar) -> None:
        """Construct a vector.

        components: Each component of the vector.
        """
        self._components = list(map(Scalar, components))

    @classmethod
    def from_scalar(cls, scalar: Scalar, dimension: int) -> VectorN:
        """Construct a vector from a scalar component.

        scalar: The single scalar component.
        dimension: The number of dimensions of this vector.

        returns: A vector with ``dimension`` components, each ``scalar``.
        """
        return cls(*([scalar] * dimension))

    @classmethod
    def basis(cls, idx: int, dimension: int) -> VectorN:
        """Construct an n-dimensional basis vector.

        idx: Basis index, e.g. 1 for y-hat.
        dimension: The number of dimensions of this vector.

        returns: A vector with 1 in the ``idx`` component and zero elsewhere.
        """
        return cls(*(int(i == idx) for i in range(dimension)))

    def __getitem__(self, idx: Union[int, slice]) -> Union[Scalar, List[Scalar]]:
        """Get numbered components."""
        return self.components[idx]

    def __setitem__(self, idx: Union[int, slice], value):
        """Set numbered components."""
        if isinstance(idx, int):
            try:
                value = Scalar(value)
            except (TypeError, ValueError):
                raise TypeError('must assign numerical value '
                                'to component') from None
        self.components[idx] = value

    def __getattr__(self, name: str) -> Union[Scalar, Tuple[Scalar]]:
        """Get named components."""
        if not all(c in NAMES for c in name):
            raise AttributeError
        try:
            if len(name) == 1:
                return self[NAMES.index(name)]
            return tuple(self[NAMES.index(c)] for c in name)
        except IndexError:
            msg = f'this vector only has {self.dimension} components'
            raise TypeError(msg) from None

    def __setattr__(self, name: str, value) -> None:
        """Set and swizzle named components."""
        if not all(c in NAMES for c in name):
            super().__setattr__(name, value)
            return
        if len(name) == 1:
            try:
                value = Scalar(value)
            except (TypeError, ValueError):
                raise TypeError('must assign numerical value '
                                'to component') from None
            try:
                self[NAMES.index(name)] = value
            except IndexError:
                msg = f'this vector only has {self.dimension} components'
                raise TypeError(msg) from None
        else:
            try:
                values = list(value)
            except TypeError:
                raise TypeError('must assign iterable to swizzled components')
            if len(name) != len(values):
                raise TypeError('must assign n values to n components')
            idxs = [NAMES.index(c) for c in name]
            for i, v in zip(idxs, values):
                self[i] = v

    def __repr__(self) -> str:
        return 'Vector2(%s)' % ', '.join(map(repr, self.components))

    __str__ = __repr__

    def __eq__(self, other) -> bool:
        if isinstance(other, VectorN):
            return self.components == other.components
        # if this is a zero vector, it is equal to any other zero
        return abs(other) == 0 == abs(self)

    # binary operators

    def __add__(self, other: SV) -> SV:
        """u.__add__(v) <=> u + v"""
        if not isinstance(other, _SV):
            return NotImplemented
        if abs(self) == 0:
            return other
        if not isinstance(other, VectorN):
            return MultiVector(self, other)
        if self.dimension != other.dimension:
            raise TypeError('cannot add vectors with different dimensions')
        return VectorN(*(a + b for a, b in zip(
            self.components, other.components)))

    def __sub__(self, other: SV) -> SV:
        """u.__sub__(v) <=> u - v"""
        return self + (-other)

    def __mul__(self, other: Union[float, VectorN]) \
            -> Union[VectorN, Blade, MultiVector]:
        """u.__mul__(v) <=> u * v <=> uv"""
        if isinstance(other, float):
            return VectorN(*(i * other for i in self.components))
        if isinstance(other, VectorN):
            return (MultiVector.from_vector(self)
                    * MultiVector.from_vector(other))
        return NotImplemented

    def __truediv__(self, other: Union[float, VectorN]) -> VectorN:
        """u.__truediv__(v) <=> u / v <=> u * ~v (for vector v)"""
        if isinstance(other, float):
            return self * (1/other)
        if isinstance(other, VectorN):
            return self * ~other
        return NotImplemented

    def __rtruediv__(self, other: float) -> VectorN:
        """u.__rtruediv__(v) <=> v / u <=> v * ~u (for scalar v)"""
        if isinstance(other, float):
            return (~self) * other
        return NotImplemented

    def __matmul__(self, other: VectorN) -> VectorN:
        """u.__matmul__(v) <=> u @ v (dot product)"""
        if not isinstance(other, VectorN):
            return NotImplemented
        if self.dimension != other.dimension:
            raise TypeError(
                'cannot multiply vectors with different dimensions')
        return (self * other + other * self) * 0.5

    def __xor__(self, other: VectorN) -> Union[Blade, MultiVector]:
        """u.__xor__(v) <=> u ^ v (outer product)"""
        if not isinstance(other, VectorN):
            return NotImplemented
        return (self * other - other * self) * 0.5

    # unary operators

    def __neg__(self) -> VectorN:
        """u.__neg__() <=> -u"""
        return VectorN(*(-i for i in self.components))

    def __pos__(self) -> VectorN:
        """u.__pos__() <=> +u <=> u normalized"""
        return self / abs(self)

    def __abs__(self) -> Scalar:
        """u.__abs__() <=> abs(u) <=> sqrt(u @ u)"""
        return (self @ self) ** .5

    def __invert__(self) -> VectorN:
        """u.__invert__() <=> ~u <=> 1/u"""
        return self / (self @ self)

class Blade(Simple):
    """u * v = Blade.from_2_vectors(u, v)"""
    bases: List[int]
    scalar: Scalar
    dimension: int

    def __new__(cls, *bases: int, dimension: int, scalar: Scalar = 1.0) -> Simple:
        bases, scalar = cls.condense_bases(bases)
        if len(bases) == 0:
            return scalar
        return super().__new__(cls)

    def __init__(self, *bases: int, dimension: int, scalar: Scalar = 1.0) -> None:
        self.dimension = int(dimension)
        self.bases, self.scalar = self.condense_bases(bases, scalar)

    @staticmethod
    def condense_bases(bases: Tuple[int], scalar: Scalar = 1.0) \
            -> Tuple[List[int], Scalar]:
        bases = list(bases)
        scalar = Scalar(scalar)
        if count_swaps(bases) % 2 == 1:
            scalar *= -1
        bases.sort()
        bases = sum(([basis] * (bases.count(basis) % 2)
                     for basis in set(bases)), [])
        return bases, scalar

    @classmethod
    def from_2_vectors(cls, u: VectorN, v: VectorN, scalar: Scalar = 1.0) \
            -> Union[Blade, MultiVector]:
        if u.dimension != v.dimension:
            raise TypeError('cannot take outer product of '
                            'vectors with different dimension')
        pairs = []
        for i, j in combinations(range(u.dimension), 2):
            pairs.append(cls(i, j,
                             scalar=u[i] * v[j] - u[j] * v[i],
                             dimension=u.dimension))
        return MultiVector(*pairs) * scalar

    def __repr__(self) -> str:
        if max(self.bases) > 3: # 0 = x, 3 = w
            bases = ' * '.join(f'e{i}' for i in self.bases)
        else:
            bases = ' * '.join(NAMES[i] for i in self.bases)
        return f'({self.scalar} * {bases})'

    __str__ = __repr__

    def __add__(self, other: SV) -> Union[Blade, MultiVector]:
        if isinstance(other, Blade):
            if self.dimension != other.dimension:
                raise TypeError(
                    'cannot combine blades of different dimensions')
            if self.bases == other.bases:
                return Blade(*self.bases, dimension=self.dimension,
                             scalar=self.scalar + other.scalar)
        return MultiVector(self, other)

    def __sub__(self, other: SV) -> Union[Blade, MultiVector]:
        return self + (-other)

    def __mul__(self, other: Union[VectorN, Simple]) -> Union[Blade, MultiVector]:
        if isinstance(other, float):
            return Blade(*self.bases, dimension=self.dimension,
                         scalar=self.scalar * other)
        if self.dimension != other.dimension:
            raise TypeError('cannot combine blades of different dimensions')
        if isinstance(other, Blade):
            return Blade(*self.bases, *other.bases,
                         dimension=self.dimension,
                         scalar=self.scalar * other.scalar)
        if isinstance(other, VectorN):
            return self * MultiVector.from_vector(other)
        return NotImplemented

    def __rmul__(self, other: Union[VectorN, Scalar]) -> Union[Blade, MultiVector]:
        if isinstance(other, float):
            return self * other
        if isinstance(other, VectorN):
            return MultiVector.from_vector(other) * self
        return NotImplemented

    def __matmul__(self, other: SV) -> Union[Blade, MultiVector]:
        """u.__matmul__(v) <=> u @ v (inner product)"""
        return (self * other + other * self) * 0.5

    def __xor__(self, other: SV) -> Union[Blade, MultiVector]:
        """u.__xor__(v) <=> u ^ v (outer product)"""
        return (self * other - other * self) * 0.5

    def __neg__(self) -> Blade:
        return self * -1

    def __abs__(self) -> Scalar:
        return self.scalar

class MultiVector:
    """Anything not simple is a multivector."""
    terms: List[Simple]

    def __new__(cls, *terms: SV) -> SV:
        terms2 = cls.condense_terms(terms)
        if len(terms2) == 1:
            return terms2[0]
        return super().__new__(cls)

    def __init__(self, *terms: SV) -> None:
        self.terms = self.condense_terms(terms)

    @staticmethod
    def condense_terms(terms: Tuple[SV]) -> Tuple[SV]:
        _terms: List[Simple] = []
        for term in terms:
            if abs(term) == 0:
                continue
            if isinstance(term, VectorN):
                term = MultiVector.from_vector(term)
            if isinstance(term, MultiVector):
                _terms.extend(term.terms)
            else:
                _terms.append(term)
        termdict: Dict[int, Dict[Tuple[int], List[Simple]]] = {}
        for term in _terms:
            termdict.setdefault(term.grade, {}).setdefault(
                tuple(term.bases), []).append(term)
        termdict2: Dict[int, Dict[Tuple[int], Simple]] = {}
        for grade, value in termdict.items():
            termdict2.setdefault(grade, {})
            for bases, value2 in value.items():
                termdict2[grade][bases] = sum(value2)
        terms2 = [value2 for value in termdict2.values()
                  for value2 in value.values() if abs(value2) != 0]
        terms2.sort(key=lambda t: (t.grade, t.bases))
        return tuple(terms2)

    @classmethod
    def from_vector(self, v: VectorN) -> MultiVector:
        """Convert a vector into a sum of basis 1-blades"""
        return sum(Blade(i, dimension=v.dimension, scalar=j)
                   for i, j in enumerate(v.components))

    def __repr__(self) -> str:
        return '(' + ' + '.join(map(repr, self.terms)) + ')'

    def __add__(self, other: SV) -> MultiVector:
        if isinstance(other, MultiVector):
            return MultiVector(*self.terms, *other.terms)
        if isinstance(other, Simple):
            return MultiVector(*self.terms, other)
        if isinstance(other, VectorN):
            return self + MultiVector.from_vector(other)
        return NotImplemented

    def __sub__(self, other: SV) -> MultiVector:
        return self + (-other)

    def __mul__(self, other: SV) -> MultiVector:
        if isinstance(other, MultiVector):
            return MultiVector(*(
                a * b for a in self.terms for b in other.terms))
        if isinstance(other, Simple):
            return MultiVector(*(
                t * other for t in self.terms))
        if isinstance(other, VectorN):
            return self * MultiVector.from_vector(other)
        return NotImplemented

    def __neg__(self) -> MultiVector:
        return self * -1

    def __abs__(self) -> Scalar:
        return (self * self) ** .5

# any number, with or without orientation
_SV = Simple, VectorN, MultiVector
SV = Union[Simple, VectorN, MultiVector]
