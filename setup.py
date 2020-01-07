#!/usr/bin/env python

import re
import setuptools

version = "0.0.5"
# with open('magition/__init__.py', 'r') as fd:
#     version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
#                         fd.read(), re.MULTILINE).group(1)

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="blxo",
    version=version,
    author="darwinqii",
    author_email="peng.qi@usask.ca",
    description="Some functions for Bent Laue Monochromators",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/darwinqii",
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'scipy>=1.2.1',
        'pathlib'
    ],entry_points={'console_scripts':['mc=blxo:mc','math_physics=blxo:math_physics']}
)