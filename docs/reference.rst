MultiVectors Module Reference
=============================

.. default-domain:: py
.. currentmodule:: multivectors

Module Attributes
-----------------

.. data:: _

    The scalar basis multivector (a 0-blade). Usage:

    >>> from multivectors import _
    >>> 2 * _
    (2.0)

.. data:: x
.. data:: y
.. data:: z
.. data:: w

    Basis vectors for the first four dimensions. Additional bases can be swizzled on the module:

    >>> from multivectors import x, y, z, xyz
    >>> x * y * z == xyz
    True

The MultiVector class
---------------------

.. autoclass:: MultiVector
    :members:
    :special-members:
    :exclude-members: __weakref__
    :member-order: bysource

Helper Functions
----------------

These functions are not really meant for exporting, but are included for completeness.

.. autofunction:: merge
.. autofunction:: count_swaps
.. autofunction:: names_to_idxs
.. autofunction:: idxs_to_idxs
.. autofunction:: idxs_to_names
.. autofunction:: condense_bases
