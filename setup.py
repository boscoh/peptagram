#!/usr/bin/env python
from setuptools import setup

description = \
"""Peptagram is a proteomics visualizer

Docs at http://github.com/boscoh/peptagram.
"""

setup(
    name='peptagram',
    version='0.1',
    author='Bosco Ho',
    author_email='boscoh@gmail.com',
    url='http://github.com/boscoh/peptagram',
    description='proteomics visualizer',
    long_description=description,
    license='BSD',
    include_package_data = True,
    install_requires=[
        'uniprot',
        'pymzml',
    ],
    packages=['peptagram',],
    scripts=[],
)