#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function
import glob
import os
import sys
import webbrowser

import tkform

import peptagram.prophet
import peptagram.proteins
from peptagram import parse






test_params = {
  'prot_xml': 'example_data/prophet/hca-lysate-16.prot.xml',
  'files_and_labels': [('example_data/prophet/hca-lysate-16.pep.xml', '16')],
  'fasta': 'example_data/prophet/HUMAN.fasta',
  'protein_error': '0.01',
  'peptide_error': '0.05',
  'great_expect': 1E-8,
  'cutoff_expect': 1E-2,
  'exclude_seqids': '',
  'include_seqids': '',
  'title': 'TPP Prophet Example Peptagram',
  'out_dir': 'peptagram-prophet',
  'include_msms': 1,
  'match_filter': 3,
}



def convert_prophet_peptagram(params, print_fn=sys.stdout.write):
  if len(params['files_and_labels']) == 0:
    raise ValueError('No files were selected.')
  parse.check_fnames(params['prot_xml'], params['fasta'])

  proteins = {}
  prot_xml = params['prot_xml'] 
  peptide_error = float(params['peptide_error'])
  protein_error = float(params['protein_error'])
  great_expect = float(params['great_expect'])
  cutoff_expect = float(params['cutoff_expect'])
  pep_xmls = [e[0] for e in params['files_and_labels']]
  size = parse.size_str(*pep_xmls)
  print_fn('Processing %s (%s)...\n' % (pep_xmls, size))
  proteins, labels = peptagram.prophet.get_proteins_and_sources(
      prot_xml, 
      pep_xmls, 
      peptide_error=peptide_error,
      protein_error=protein_error,
      good_expect=great_expect,
      cutoff_expect=cutoff_expect)
  labels = map(parse.basename, labels)
  peptagram.proteins.filter_proteins(proteins, params)

  peptagram.proteins.make_graphical_comparison_visualisation({
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': ['expect=%s' % great_expect, 'expect=%s' % cutoff_expect],
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

    self.title('Create TPP Prophet Peptagram')

    self.push_text("Create TPP Prophet Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED INPUT DATA", 16)
    self.push_labeled_param(
        'fasta', 'Protein sequences in fasta', 'sequences.fasta', load_file_text='select .fasta')
    self.push_labeled_param(
        'prot_xml', 'ProteinProphet .prot.xml', 'job.prot.xml', load_file_text='select .prot.xml')
    self.push_text("PeptideProphet .pep.xml files; drag arrow to reorder; edit labels for peptagram")
    self.push_file_list_param(
        'files_and_labels', '+ .pep.xml files', is_label=False)

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("PEPTAGRAM PARAMETERS", 16)
    self.push_labeled_param(
        'title', 'Title', 'Prophet Peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-prophet', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL PARAMETERS", 16)

    self.push_labeled_param('peptide_error', 'PeptideProphet false-positive-error cutoff', '0.05')
    self.push_labeled_param('protein_error', 'ProteinProphet false-positive-error cutoff', '0.01')
    self.push_labeled_param('great_expect', 'Good PSM expect', '1E-8')
    self.push_labeled_param('cutoff_expect', 'Cutoff PSM expect', '1E-2')

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
    index = convert_prophet_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)



if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_prophet_peptagram(test_params)
  else:
    form = PeptagramForm(800, -50)
    form.mainloop()



