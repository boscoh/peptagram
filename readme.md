

# PEPTAGRAM

`peptagram` generates a single-page HTML5 web-app to visualize proteomics analyses for a graphical comparison between multiple proteomics experiment ([example][example]).

[example]:http://boscoh.github.io/peptagram/examples/multiple/index.html

The resulting HTML5 web-app is:

- instantly available: *runs on any modern webbrowser*
- cross-platform: *everyone has a webbrowser*
- self-contained: *just zip directory and send by email*
- publishable: *installable on any server, even on Dropbox*

More information at <http://boscoh.github.io/peptagram>.


## Installation

`peptagram` consists of a set of python scripts that converts proteomics data into an HTML5 visualisation. 

`peptagram` requries python 2.7. If your system does not have it, you should download [python](https://www.python.org/downloads/) and the python package installer [pip](https://pip.pypa.io/en/latest/installing.html). 

Download and unzip the package from <https://github.com/boscoh/peptagram/archive/master.zip>.

**For Developers**. To use `peptagram` as a module, you should:

    > pip install peptagram

Or download the package and run:

    > python setup.py install

`peptagram` has two python dependencies: 

  1. [pymzml](https://github.com/pymzml/pymzML) to read `.mzML` files
  2. [tkform](https://github.com/boscoh/tkform) to generate the GUI

To install them:

    > pip install pymzml tkform



## Usage

1. User guide to the HTML5 visualisation <http://boscoh.github.io/peptagram/vizhelp.html>.
2. Examples for generating visualisations from standard formats <http://boscoh.github.io/peptagram/vizgen.html>.
3. Programming API to use custom data <http://boscoh.github.io/peptagram/api.html>.
4. Source code and open-source contributions <http://github.com/boscoh/peptagram>.


## Credits

(c) 2013, 2015. Bosco K. Ho.  
Developed by Bosco Ho at the Monash Proteomics Facility.  
Original concept by Oded Kleifeld.  
Contributions by Rob Goode.  


## Changelog

0.5 (June 2015)

- built GUIs for all the scripts
- more responsive webapps
- can combine peptagrams
- streamlined parsers and support for more MS search types

0.2.1 (2013)

- separated the delete_empty_proteins out in the proteins module
