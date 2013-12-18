#!/usr/bin/env python
from setuptools import setup

setup(
    name='peptagram',
    version='0.1',
    author='Bosco Ho',
    author_email='boscoh@gmail.com',
    url='http://boscoh.github.io/peptagram',
    description='generates HTML5 visualization for proteomics analyses',
    long_description="Docs at http://boscoh.github.io/peptagram",
    license='BSD',
    include_package_data = True,
    install_requires=[
        'uniprot',
        'pymzml',
    ],
    packages=['peptagram',],
    scripts=[],
)