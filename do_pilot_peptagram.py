#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-


from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform


import peptagram.proteins
import peptagram.pilot
from peptagram import parse

import tkform


test_params = {
  'files_and_labels': [('example_data/pilot/DPrep1_pilot4.txt', 'summary')],
  'fasta': 'example_data/pilot/HUMAN.fasta',
  'exclude_seqids': '',
  'include_seqids': '',
  'title': 'Protein Pilot Example Petagram',
  'out_dir': 'peptagram-pilot',
  'include_msms': 0,
  'match_filter': 3,
}


def convert_pilot_to_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')

  size = parse.size_str(params['fasta'])
  print_fn("Using sequences from %s (%s)...\n" % (params['fasta'], size))

  proteins = {}
  labels = []
  for fname, label in params['files_and_labels']:
    labels.append(label)
    size = parse.size_str(fname)
    print_fn("Processing %s (%s)...\n" % (fname, size))
    these_proteins = peptagram.pilot.get_proteins(fname)
    proteins = peptagram.proteins.merge_two_proteins(
        proteins, these_proteins)

  peptagram.proteins.filter_proteins(proteins, params)

  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': [0, 1],
      'mask_labels': [],
      'out_dir': params['out_dir'],
  })

  html = os.path.join(params['out_dir'], 'index.html')
  size = parse.size_str(params['out_dir'])
  print_fn('Successfully built peptagram (%s): %s\n' % (size, html))
  return os.path.abspath(html)


class PeptagramForm(tkform.Form):
  """
  Application window for Tkinter.
  """ 
  def __init__(self, width=700, height=800, parent=None):
    tkform.Form.__init__(self, parent, width, height)
    self.title('Create ProteinPilot Peptagram')

    self.push_text("Create ProteinPilot Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED INPUT DATA", 16)
    self.push_labeled_param(
        'fasta', 'Load sequences for proteins', 'fasta', load_file_text='select .fasta')
    self.push_text("Load ProteinPilot peptide summaries; drag arrow to reorder; edit labels for peptagram")
    self.push_file_list_param(
        'files_and_labels', '+ peptide summary .csv/.txt file')
    
    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("PEPTAGRAM PARAMETERS", 16)
    self.push_labeled_param(
        'title', 'Peptagram Title', 'Protein Pilot peptagram')
    self.push_labeled_param(
        'out_dir', 'Output Directory', 'peptagram-pilot', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL FILTERS", 16)

    self.push_labeled_param(
        'exclude_seqids', 'Filename of seqids to exclude',
         load_file_text='select')
    self.push_labeled_param(
        'include_seqids', 'Filename of seqids to include',
         load_file_text='select')

    self.push_text("Include matches:")
    self.push_radio_param(
        'match_filter',
        ['All matches', 
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
    index = convert_pilot_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)


if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_pilot_to_peptagram(test_params)
  else:
    form = PeptagramForm(800, -150)
    form.mainloop()



