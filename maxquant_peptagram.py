#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform

import peptagram.proteins
import peptagram.maxquant
from peptagram import parse


test_params = {
  'files_and_labels': [('example_data/maxquant/txt', 'txt')],
  'fasta': 'example_data/maxquant/yeast_orf_trans_all_05-Jan-2010.fasta',
  'exclude_seqids': '',
  'include_seqids': '',
  'great_expect': '1E-8',
  'cutoff_expect': '1E-2',
  'title': 'Maxquant Example Peptagram',
  'out_dir': 'peptagram-maxquant',
  'include_msms': 1,
  'match_filter': 3,
}


def convert_maxquant_to_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')

  parse.check_fnames(params['fasta'])

  great_expect = float(params['great_expect'])
  cutoff_expect = float(params['cutoff_expect'])

  proteins = {}
  labels = []
  for fname, label in params['files_and_labels']:
    size = parse.size_str(fname)
    print_fn("Processing %s (%s)...\n" % (fname, size))
    labels.extend(label)
    these_proteins, sources = \
        peptagram.maxquant.get_proteins_and_sources(
            fname, great_expect=great_expect, cutoff_expect=cutoff_expect)
    proteins = peptagram.proteins.merge_two_proteins(
        proteins, these_proteins)

  peptagram.proteins.filter_proteins(proteins, params)

  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': ['PEP=%s'%params['great_expect'], 
                      'PEP=%s'%params['cutoff_expect']],
      'mask_labels': [],
      'out_dir': params['out_dir'],
  })

  html = os.path.join(params['out_dir'], 'index.html')
  size = parse.size_str(html)
  print_fn('Successfully built peptagram (%s): %s\n' % (size, html))
  return os.path.abspath(html)




class PeptagramForm(tkform.Form):
  """
  Application window for Tkinter.
  """ 
  def __init__(self, width=700, height=800, parent=None):
    tkform.Form.__init__(self, parent, width, height)
    self.title('Create MaxQuant Peptagram')

    self.push_text("Create MaxQuant Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED PARAMTERS", 20)
    self.push_text("Maxquant summary txt directory; drag arrow to reorder; edit labels for peptagram")

    self.push_dir_list_param(
        'files_and_labels', '+ txt directory')
    self.push_labeled_param(
        'fasta', 'Protein sequences', 'fasta', load_file_text='select .fasta')
    self.push_labeled_param(
        'title', 'Peptagram title', 'Maxquant peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-maxquant', load_dir_text='select directory')

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
    self.push_labeled_param('great_expect', 'Good PEP', '1E-8')
    self.push_labeled_param('cutoff_expect', 'Cutoff PEP', '1E-2')
    self.push_checkbox_param('include_msms', 'Include MS/MS')
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
    index = convert_maxquant_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)



if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_maxquant_to_peptagram(test_params)
  else:
    form = PeptagramForm(800, -50)
    form.mainloop()



