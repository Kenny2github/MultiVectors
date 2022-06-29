import re
from setuptools import setup

with open('multivectors.py') as f:
    text = f.read()
match = re.search('"""(.+?)"""', text, re.S)
if match is None:
    raise ValueError('Could not find module docstring in multivectors.py')
longdesc = match.group(1).strip()
with open('README.md', 'w') as f:
    f.write(longdesc)

setup()
