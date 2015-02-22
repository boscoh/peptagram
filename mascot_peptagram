#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform

import peptagram.proteins
from peptagram import parse
import peptagram.mascot





test_params = {
  'files_and_labels': [('example_data/mascot/F022043.dat', 'F022043')],
  'exclude_seqids': '',
  'include_seqids': '',
  'great_ionscore': 80,
  'cutoff_ionscore': 0,
  'title': 'Mascot Peptagram Example',
  'out_dir': 'example_data/peptagram-mascot',
  'fasta': 'example_data/mascot/HUMAN.fasta',
  'include_msms': 1,
  'match_filter': 3,
}


def convert_mascot_to_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')
  great_ionscore = float(params['great_ionscore'])
  cutoff_ionscore = float(params['cutoff_ionscore'])
  proteins = {}
  labels = []
  for fname, label in params['files_and_labels']:
    labels.append(label)
    size = parse.size_str(fname)
    print_fn("Processing %s (%s)...\n" % (fname, size))
    these_proteins = peptagram.mascot.get_proteins(
        fname, great_ionscore, cutoff_ionscore)
    proteins = peptagram.proteins.merge_two_proteins(
        proteins, these_proteins)
  peptagram.proteins.filter_proteins(proteins, params)
  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': [great_ionscore, cutoff_ionscore],
      'mask_labels': [],
      'out_dir': params['out_dir'],
  })
  html = os.path.join(params['out_dir'], 'index.html')
  size = parse.size_str(html)
  print_fn('Successfully built peptagram (%s): %s\n' % (size, html))
  return os.path.abspath(html)


class PeptagramForm(tkform.Form):

  def __init__(self, width=700, height=800, parent=None):
    tkform.Form.__init__(self, parent, width, height)

    self.title('Create Mascot Peptagram')

    self.push_text("Create Mascot Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED PARAMETERS", 20)
    self.push_text("Mascot .dat files; drag arrow to reorder; edit labels for peptagram")
    self.push_file_list_param('files_and_labels', '+ .dat files')
    self.push_labeled_param(
        'fasta', 'Protein sequences', 'sequences.fasta', load_file_text='select .fasta')
    self.push_labeled_param(
        'title', 'Peptagram title', 'Mascot Peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-mascot', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL FILTERS", 20)

    self.push_labeled_param(
        'exclude_seqids', 'Text file of excluded seqids',
         load_file_text='select')
    self.push_labeled_param(
        'include_seqids', 'Text file of included seqids',
         load_file_text='select')
    self.push_labeled_param('great_ionscore', 'Good ionscore', '80')
    self.push_labeled_param('cutoff_ionscore', 'Cutoff ionscore', '0')
    self.push_checkbox_param('include_msms', 'Include MS/MS')
    self.push_labeled_param(
        'n_peak', 'Number of peaks in spectrum', '50')
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

    self.push_text("OUTPUT", 20)
    self.push_submit()
    self.push_output()

  def run(self, params):
    index = convert_mascot_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)


if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_mascot_to_peptagram(test_params)
    sys.exit(1)
  else:
    form = PeptagramForm(800, -50)
    form.mainloop()




