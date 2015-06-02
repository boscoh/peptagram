#!/usr/bin/env python

"""
Assuming prototyping is done in the peptagram here,
copies the modified js back into the master python
template directory.
"""

import os
import shutil
import peptagram


fname_pairs = """
copyback_comparison.py templates/comparison/
index.html  templates/comparison/
js/peptagram.js  templates/js/
js/canvaswidget.js  templates/js/
js/datacontroller.js  templates/js/
js/mass.js  templates/js/
js/matchlist.js  templates/js/
js/proteinlist.js  templates/js/
js/spectrumwidget.js  templates/js/
js/util.js  templates/js/
"""


def copy_peptagram_files(in_dir, out_dir):
  fnames = fname_pairs.split()

  pairs = []
  for i in range(0, len(fnames), 2):
    in_fname, out_rel_dir = fnames[i], fnames[i+1]
    in_fname = os.path.join(in_dir, in_fname)
    out_fname = os.path.join(
        out_dir, out_rel_dir, os.path.basename(in_fname))
    pairs.append((in_fname, out_fname))

  max_len = max([len(i) for i,o in pairs])
  for i, o in pairs:
    format_str = '%% %ds -> %%s' % max_len
    print str(i).ljust(max_len) + ' -> ' + o
    # assert os.path.isfile(i)
    # assert os.path.isfile(o)
    shutil.copy(i, o)


copy_peptagram_files('.', os.path.dirname(peptagram.__file__))







