#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform

import peptagram.xtandem
import peptagram.proteins
from peptagram import parse


test_params = {
  'files_and_labels': [('example_data/xtandem/Seq23282_E1O1.tandem', 'Seq23282')],
  'exclude_seqids': '',
  'include_seqids': '',
  'title': 'X!Tandem Example Peptagram',
  'out_dir': 'peptagram-xtandem',
  'fasta': '',
  'great_expect': 1E-8,
  'okay_expect': 1E-4,
  'cutoff_expect': 1E-2,
  'include_msms': 1,
  'match_filter': 3,
  'n_peak': 50,
}


def convert_xtandem_to_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')

  n_peak = int(params['n_peak'])
  great_expect = float(params['great_expect'])
  cutoff_expect = float(params['cutoff_expect'])

  proteins = {}
  labels = []
  for fname, label in params['files_and_labels']:
    labels.append(label)
    size = parse.size_str(fname)
    print_fn("Processing %s (%s)...\n" % (fname, size))
    these_proteins = peptagram.xtandem.get_proteins(
        fname, 
        n_peak=n_peak,
        good_expect=great_expect,
        cutoff_expect=cutoff_expect)
    proteins = peptagram.proteins.merge_two_proteins(
        proteins, these_proteins)

  peptagram.proteins.filter_proteins(proteins, params)

  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': [great_expect, cutoff_expect],
      'mask_labels': [],
      'out_dir': params['out_dir'],
  })

  html = os.path.join(params['out_dir'], 'index.html')
  size = parse.size_str(params['out_dir'])
  print_fn('Successfully built peptagram (%s): %s\n' % (size, html))
  return os.path.abspath(html)



class PeptagramForm(tkform.Form):

  def __init__(self, width=700, height=800, parent=None):
    tkform.Form.__init__(self, parent, width, height)

    self.title('Create X!Tandem Peptagram')

    self.push_text("Create X!Tandem Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED INPUT DATA", 16)
    self.push_text(u"Load X!Tandem files; drag \u2630 to reorder; edit labels")
    self.push_file_list_param(
        'files_and_labels', '+ .tandem files')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("PEPTAGRAM PARAMETERS", 16)
    self.push_labeled_param(
        'title', 'Title', 'X!Tandem Peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-xtandem', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL FILTERS", 16)

    self.push_labeled_param(
        'exclude_seqids', 'Text file of excluded seqids',
         load_file_text='select')
    self.push_labeled_param(
        'include_seqids', 'Text file of included seqids',
         load_file_text='select')

    self.push_labeled_param('great_expect', 'Good expect', '1E-8')
    self.push_labeled_param('cutoff_expect', 'Cutoff expect', '1E-2')

    self.push_checkbox_param('Include_msms', 'include MS/MS')
    self.push_labeled_param(
        'n_peak', 'Number of peaks for spectrum', '50')
    self.push_text("Include matches:")
    self.push_radio_param(
        'match_filter',
        ['All peptides', 
         'Tryptic peptides', 
         'Semi-tryptic peptides', 
         'Modified peptides'])

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("RESULTS", 16)
    self.push_submit()
    self.push_output()


  def run(self, params):
    index = convert_xtandem_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)



if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_xtandem_to_peptagram(test_params)
  else:
    form = PeptagramForm(800, -150)
    form.mainloop()



