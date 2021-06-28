import re
from setuptools import setup

with open('multivectors.py') as f:
    text = f.read()
longdesc = re.search('"""(.+?)"""', text, re.S).group(1).strip()
version = re.search(r"__version__\s*=\s*'([^']+)'", text).group(1).strip()
with open('README.md', 'w') as f:
    f.write(longdesc)

setup(
    name='multivectors',
    version=version,
    license='Apache',
    description='Compute blades and multivectors in arbitrary-dimensional space.',
    long_description=longdesc,
    long_description_content_type='text/markdown',
    author='Ken Hilton',
    url='https://github.com/Kenny2github/MultiVectors',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: Apache Software License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Typing :: Typed'
    ],
    keywords='multivector geometric algebra',
    py_modules=['multivectors'],
    python_requires='>=3.7'
)
