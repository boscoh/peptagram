

# PEPTAGRAM 


`peptagram` generates a single-page HTML5 web-app to visualize proteomics analyses for:

1. an overview of a proteomics experiment ([example][example1]),
2. a comparison between multiple proteomics experiment ([example][example2]).

[example1]:http://boscoh.github.io/peptagram/examples/overview/index.html
[example2]:http://boscoh.github.io/peptagram/examples/comparison/index.html

The resulting HTML5 web-app is:

- instantly available: *runs on any modern webbrowser*
- cross-platform: *everyone has a webbrowser*
- self-contained: *just zip directory and send by email*
- publishable: *installable on any server, even on Dropbox*

More information at <http://boscoh.github.io/peptagram>.


## Installation

`peptagram` consists of a set of python scripts that converts proteomics data into an HTML5 visualisation. 

To generate the visualisations, you must have python installed. If you have the python installer pip, then:

    > pip install peptagram

Otherwise, download and unzip the package from <https://github.com/boscoh/peptagram/archive/master.zip>

And install from the package directory:

    > python setup.py install

`peptagram` has two python dependencies: 

  1. [pymzml](https://github.com/pymzml/pymzML) to read `.mzML` files
  2. [uniprot](https://github.com/boscoh/uniprot) to get protein sequences from <http://uniprot.org>. 

These should be installed automatically from the above scripts, but if that fails, you may need to install them manually.


## Usage

1. User guide to the HTML5 visualisation <http://boscoh.github.io/peptagram/vizhelp.html>.

2. Examples for generating visualisations from standard xformats <http://boscoh.github.io/peptagram/vizgen.html>.

3. Programming API to use custom data <http://boscoh.github.io/peptagram/api.html>.

4. Source code and open-source contributions <http://github.com/boscoh/peptagram>.



##

Copyright (c) 2013, Bosco K. Ho. 
Supported by Monash University.
