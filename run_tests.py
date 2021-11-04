import sys
import doctest
import multivectors

fails = doctest.testmod(multivectors)[0]
fails += doctest.testfile('README.md')[0]
fails += doctest.testfile('docs.md')[0]
if fails > 0:
    sys.exit(1)