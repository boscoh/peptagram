#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform

import peptagram.morpheus
import peptagram.mzml
import peptagram.proteins
from peptagram import parse



test_params = {
  'files_and_labels': [('example_data/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv', 'KO1')],
  'mzmls_and_labels': [('example_data/morpheus/OK20130822_MPProtomap_KO1.mzML', '')],
  'modifications': 'example_data/morpheus/modifications.tsv',
  'exclude_seqids': '',
  'include_seqids': '',
  'title': 'Morpheus Example Peptagram',
  'out_dir': 'peptagram-morpheus',
  'include_msms': 1,
  'match_filter': 3,
  'n_peak': 50,
  'q_cutoff': 10,
  'q_good': 0,
}



def convert_morpheus_to_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')

  q_good = float(params['q_good'])
  q_cutoff = float(params['q_cutoff'])

  modifications = params['modifications']
  proteins = {}
  labels = []
  for entry in params['files_and_labels']:
    fname = entry[0]
    size = parse.size_str(fname)
    print_fn("Processing %s (%s)...\n" % (fname, size))
    protein_group = fname.replace('PSMs', 'protein_groups')
    print_fn("Inferring protein group: %s\n" % protein_group)
    these_proteins, these_sources = \
        peptagram.morpheus.get_proteins_and_sources(
            protein_group, fname, modifications, q_good, q_cutoff)
    proteins = peptagram.proteins.merge_two_proteins(
        proteins, these_proteins)
    labels.extend(map(parse.basename, these_sources))

  n_peak = int(params['n_peak'])
  for i_source, (mzml, label) in enumerate(params['mzmls_and_labels']):
    peptagram.mzml.load_mzml(
        proteins, i_source, mzml, n_peak)

  peptagram.proteins.filter_proteins(proteins, params)

  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': ['Q=%s' % q_good, 'Q=%s' % q_cutoff],
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

    self.title('Create Morpehus Peptagram')

    self.push_text("Create Morpehus Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED INPUT DATA", 16)
    self.push_text("Morpheus .PSMs.tsv files; drag arrow to reorder; edit labels for peptagram")
    self.push_file_list_param(
        'files_and_labels', '+ .PSMs.tsv files', is_label=False)


    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("PEPTAGRAM PARAMETERS", 16)
    self.push_labeled_param(
        'title', 'Title', 'Morpheus Peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-morpheus', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL PARAMETERS TO DISPLAY SPECTRA", 16)
    self.push_text("The .mzML file that produced the search")
    self.push_file_list_param('mzmls_and_labels', '+ .mzML')
    self.push_labeled_param(
        'n_peak', 'Number of peaks for spectrum', '50')
    self.push_text('To show matched peaks of modified peptides, select "Modifications.tsv":')
    self.push_labeled_param(
        'modifications', '', '', load_file_text='select file')
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
    self.push_labeled_param('q_good', 'Good Q', '0')
    self.push_labeled_param('q_cutoff', 'Cutoff Q', '10')
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
    index = convert_morpheus_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)



if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_morpheus_to_peptagram(test_params)
  else:
    form = PeptagramForm(800, -50)
    form.mainloop()




