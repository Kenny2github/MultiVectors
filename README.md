# multivectors
Compute blades and multivectors in arbitrary-dimensional space.

## Installation
```
pip install multivectors
```

## Usage
```py
import math
from multivectors import x, y, z

v = 2*x + 3*y + 4*z
print(v.rotate(math.pi/2, x * y))
```
Output:
```
(-3.00x + 2.00y + 4.00z)
```