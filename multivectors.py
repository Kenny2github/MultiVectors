from __future__ import annotations
from itertools import combinations
import math
from numbers import Real
from typing import Dict, List, Tuple, Union

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

def names_to_idxs(name: str) -> List[int]:
    """Convert swizzled basis vector names into generalized basis indices.

    Examples:
    - ``names_to_idxs('xyw') == [0, 1, 3]``
    - ``names_to_idxs('e_1e2_z') == [0, 1, 2]``
    - ``names_to_idxs('_') == []``
    """
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
            while j < nl and '0' <= name[j] <= '9':
                j += 1
            try:
                idxs.append(int(name[i:j]) - 1)
            except ValueError:
                raise TypeError('cannot have empty eN notation') from None
            i = j
        else:
            i += 1
    return idxs

def idxs_to_idxs(idxs: Index) -> List[int]:
    """Convert multiple possible ways to specify multiple indices.

    Examples:
    - ``idxs_to_idxs(slice(None, 5, None)) == [0, 1, 2, 3, 4]``
    - ``idxs_to_idxs((1, 3, 4)) == [1, 3, 4]``
    - ``idxs_to_idxs(1) == [1]``
    """
    if isinstance(idxs, int):
        return [idxs]
    if isinstance(idxs, slice):
        if idxs.stop is None:
            raise TypeError('cannot have infinite bases')
        return list(range(*idxs.indices(idxs.stop)))
    return idxs

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
    - ``Blade.xy == Blade(0, 1)``
    - ``Blade.e1_e3 == Blade(0, 2)``
    - ``Blade._ == Blade()``
    And indices can be combined:
    - ``Blade[1, 3] == Blade(1, 3)``
    - ``Blade[1:3] == Blade(1, 2)``
    """
    bases: Tuple[int, ...]
    scalar: Real

    _insts = {}

    @property
    def grade(self) -> int:
        """The k of "k-blade". 1 for basis vectors, 0 for scalars."""
        return len(self.bases)

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
        If so, promote it to a 0-blade.
        Otherwise, just return the object.
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
        - ``condense_bases((1, 1, 2, 1, 2), 2.0) == ((1,), -2.0)``
        - ``condense_bases((1, 2, 1, 2), 1.5) == ((), -1.5)``
        - ``condense_bases((2, 1, 3, 2, 3, 3), 1.0) == ((1, 3), 1.0)``
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
        """
        if not self.bases:
            return str(self.scalar)
        if max(self.bases) > 3:
            r = ''.join(f'e{i}' for i in self.bases)
        else:
            r = ''.join(NAMES[i] for i in self.bases)
        return '%s%s' % (self.scalar, r)

    # Binary operators

    def __add__(self, other: MVV) -> MV:
        """Add a blade and another object.

        Returns: 0-blade if this blade has grade 0 and the object is scalar.
        Returns: k-blade with scalars summed if the object is a k-blade.
        Returns: MultiVector summing this object and the other otherwise.
        """
        if isinstance(other, Real) and self.bases == ():
            return Blade(scalar=self.scalar + other)
        if isinstance(other, Blade) and self.bases == other.bases:
            return Blade(*self.bases, scalar=self.scalar + other.scalar)
        if isinstance(other, _MVV):
            return MultiVector.from_terms(self, other)
        return NotImplemented

    def __radd__(self, other: Real) -> Blade:
        """Support adding blades on the right side of scalars."""
        if isinstance(other, Real):
            return self + other # scalar addition is commutative
        return NotImplemented

    def __sub__(self, other: MVV) -> MV:
        """Subtracting is adding the negation."""
        return self + (-other)

    def __rsub__(self, other: Real) -> Blade:
        """Support subtracting blades from scalars."""
        return other + (-self)

    def __mul__(self, other: Simple) -> Blade:
        """Multiply a blade and another object.

        Returns: k-blade scaled by object if object is scalar.
        Returns: The geometric product, if both objects are blades.
        """
        if isinstance(other, Real):
            return Blade(*self.bases, scalar=self.scalar * other)
        if isinstance(other, Blade):
            return Blade(*self.bases, *other.bases,
                         scalar=self.scalar * other.scalar)
        return NotImplemented

    def __rmul__(self, other: Real) -> Blade:
        """Support multiplying blades on the right side of scalars."""
        if isinstance(other, Real):
            return self * other # scalar multiplication is commutative
        return NotImplemented

    def __matmul__(self, other: Simple) -> Simple:
        """Get the inner (dot) product of two objects.

        Returns: The scaled blade when dotting with a scalar.
        Returns: The inner product of the two blades, which can be scalar 0.
        """
        return (self * other) % abs(self.grade - Blade.check(other).grade)

    dot = inner = __matmul__

    def __truediv__(self, other: Simple) -> Blade:
        """Divide two objects.

        Returns: The scaled blade when dividing by a scalar.
        Returns: The geometric product of this blade and the other's inverse.
        """
        if isinstance(other, Real):
            return self * (1 / other)
        if isinstance(other, Blade):
            return self * other / (abs(other) * abs(other))
        return NotImplemented

    def __rtruediv__(self, other: Real) -> Blade:
        """Divide a scalar by a blade."""
        return other * self / (abs(self) * abs(self))

    def __mod__(self, other: int) -> Simple:
        """The choose operator - returns the sum of all blades of grade k.
        When applied to a blade, returns the blade if k is the blade's grade.
        Returns scalar 0 otherwise.
        """
        return self if self.grade == other else 0.0

    choose = __mod__

    def __pow__(self, other: int, modulo: int = None) -> Simple:
        """A blade raised to an integer power.
        A ** n = A * A * A * ... * A n times.
        The optional ternary power uses the choose operator afterwards.
        """
        if not isinstance(other, int):
            if not self.bases:
                return float(self) ** other
            return NotImplemented
        result = 1
        for _ in range(other):
            result *= self
        if modulo is not None:
            if not isinstance(modulo, int):
                return NotImplemented
            return result % modulo
        return result

    def __rpow__(self, other: Real) -> Simple:
        """A real number raised to a blade power.
        x ** A = e ** (A ln x) = e ** (I * a ln x)
        = cos(a ln x) + sin(a ln x) * I
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
        Returns: The outer product of the two blades, which can be scalar 0.
        """
        return (self * other) % (self.grade + other.grade)

    wedge = outer = __xor__

    # Unary operators

    def __neg__(self) -> Blade:
        """The negation of a blade is the blade with its scalar negated."""
        return self * -1

    def __abs__(self) -> Real:
        """Get the scalar of a blade."""
        return float(self.scalar)

    magnitude = __abs__

    def __invert__(self) -> Blade:
        """A normalized blade is one with a magnitude of 1."""
        return Blade(*self.bases)

    normalize = __invert__

    def __float__(self) -> float:
        """Convert the scalar of a blade to a float."""
        if self.bases == ():
            return float(self.scalar)
        return NotImplemented

class MultiVector:
    """Two or more incompatible blades, summed together.

    The bare constructor is not meant for regular use.
    Use the factory ``MultiVector.from_terms()`` instead.

    Basis vector names can be swizzled on instances:
    - ``(x + y).x == 1.0``
    - ``(x * y + z).xy == 1.0``
    - ``(x + y).e3 == 0.0``
    And indices can be combined:
    - ``(x + y)[0] == 1.0``
    - ``(x * y + z)[0, 1] == 1.0``
    - ``(x + y)[2] == 0.0``
    """

    termdict: Dict[Tuple[int, ...], Real]

    @property
    def terms(self) -> Tuple[Blade, ...]:
        """Get a sequence of blades comprising this multivector."""
        keys = sorted(self.termdict.keys(), key=lambda b: (len(b), b))
        return tuple(Blade(*key, scalar=self.termdict[key]) for key in keys)

    def __init__(self, termdict: Dict[Tuple[int, ...], Real]):
        self.termdict = termdict.copy()

    @classmethod
    def from_terms(cls, *terms: MV) -> MV:
        """Create a multivector from a sequence of terms.
        This may return something other than a MultiVector
        if it is the only term or there are no terms.
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
        idxs = tuple(names_to_idxs(name))
        return self[idxs]

    def __getitem__(self, idxs: Index) -> Union[Real, Tuple[Real]]:
        """Support index swizzling."""
        idxs = tuple(idxs_to_idxs(idxs))
        return self.termdict.get(idxs, 0.0)

    def __repr__(self) -> str:
        """Return a representation of this multivector.
        Depending on the global namespace, this may be eval()-able.
        """
        return '(' + ' + '.join(map(repr, self.terms)) + ')'

    def __str__(self) -> str:
        """Return a representation of this multivector suited for showing."""
        return '(' + ' + '.join(map(str, self.terms)) + ')'

    # Binary operators

    def __add__(self, other: MVV) -> MV:
        """Add a multivector and another object.

        Returns: The sum of the terms of this multivector and the other.
        Returns: The other object added to this multivector's terms.
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*self.terms, *other.terms)
        if isinstance(other, _Simple):
            return self.from_terms(*self.terms, other)
        return NotImplemented

    def __radd__(self, other: MVV) -> MV:
        """Support adding multivectors on the right side of objects."""
        if isinstance(other, MultiVector):
            return self.from_terms(*other.terms, *self.terms)
        if isinstance(other, _Simple):
            return self.from_terms(other, *self.terms)
        return NotImplemented

    def __sub__(self, other: MVV) -> MV:
        """Subtracting is adding the negation."""
        return self + (-other)

    def __rsub__(self, other: MVV) -> MV:
        """Support subtracting multivectors from objects."""
        return other + (-self)

    def __mul__(self, other: MVV) -> MV:
        """Multiply a multivector and another object.

        Returns: (a+b)*(c+d)=a*c+a*d+b*c+b*d for multivectors (a+b) and (c+d)
        Returns: (a+b)*v = a*v + b*v for multivector (a+b) and simple v
        """
        if isinstance(other, MultiVector):
            return self.from_terms(*(a * b for b in other.terms
                                     for a in self.terms))
        if isinstance(other, _Simple):
            return self.from_terms(*(t * other for t in self.terms))
        return NotImplemented

    def __rmul__(self, other: Simple) -> MV:
        """Support multiplying multivectors on the right side of simples."""
        if isinstance(other, _Simple):
            return self.from_terms(*(other * t for t in self.terms))
        return NotImplemented

    def __matmul__(self, other: MVV) -> MV:
        """Get the inner (dot) product of two objects.

        Returns: (a+b)@(c+d)=a@c+a@d+b@c+b@d for multivectors (a+b) and (c+d)
        Returns: (a+b)@v = a@v + b@v for multivector (a+b) and simple v
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

    def __truediv__(self, other: Simple) -> MV:
        """Divide two objects.

        Returns: (a+b)/v = a/v + b/v for multivector (a+b) and simple v
        """
        if isinstance(other, _Simple):
            return self.from_terms(*(t / other for t in self.terms))
        return NotImplemented

    def __mod__(self, grade: int) -> MVV:
        """The choose operator - returns the sum of all blades of grade k.

        Examples:
        - (a+bx+cy+dxy) % 1 = bx+cy
        - (a+b+cx+dxy+exy) % 2 = dxy+exy
        - (a+b+cx+dxy+exy) % 0 = a+b
        """
        if not isinstance(grade, int):
            return NotImplemented
        bs = set(basis for bases in self.termdict.keys() for basis in bases)
        return sum(self[bases] for bases in combinations(bs, grade))

    choose = __mod__

    def __pow__(self, other: int, modulo: int = None) -> MVV:
        """A multivector raised to an integer power.
        V ** n = V * V * V * ... * V n times.
        The optional ternary power uses the choose operator afterwards.
        """
        if not isinstance(other, int):
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

    # Unary operators

    def __neg__(self) -> MultiVector:
        """The negation of a multivector is the negation of all its terms."""
        return self.from_terms(*(-t for t in self.terms))

    def __abs__(self) -> Real:
        """The magnitude of a multivector is the square root of
        the sum of the squares of its terms.
        """
        return math.fsum(t * t for t in self.terms) ** .5

    magnitude = __abs__

    def __invert__(self) -> MultiVector:
        """A normalized multivector is one scaled down by its magnitude."""
        return self / abs(self)

    normalize = __invert__

    # Actual methods

    def rotate(self, angle: Real, plane: Blade) -> MultiVector:
        """Rotate this multivector by angle in rads around the blade plane.

        angle: Angle to rotate by, in radians.
        plane: Blade representing basis plane to rotate through.
        """
        power = plane * angle / 2
        return (math.e ** (-power)) * self * (math.e ** power)

    def angle_to(self, other: MultiVector) -> MultiVector:
        """Get the angle between this multivector and another."""
        return math.acos((self @ other) / (abs(self) * abs(other)))

_ = Blade._
x = Blade.x
y = Blade.y
z = Blade.z
w = Blade.w

Simple = Union[Real, Blade]
_Simple = Real, Blade
MV = Union[Blade, MultiVector]
MVV = Union[Real, Blade, MultiVector]
_MVV = Real, Blade, MultiVector
Index = Union[int, Tuple[int, ...], slice]
