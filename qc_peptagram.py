#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

from __future__ import print_function

import glob
import os
import sys
import webbrowser
import csv
from collections import defaultdict

import tkform

import peptagram.proteins
from peptagram import parse





test_params = {
  'skyline_csv': 'TotalTransitionsCombinedLibImport3.csv',
  'exclude_seqids': '',
  'include_seqids': '',
  'fasta': 'test.fasta',
  'title': 'EQ Peptagram',
  'out_dir': 'example_data/peptagram-qc',
  'include_msms': 0,
  'match_filter': 0,
}


def add_source_to_proteins(proteins):
  for seqid, protein in proteins.items():
    protein['sources'].append({ 'matches':[] })


def average_same_peptides_in_source(source):
    matches = source['matches']
    matches_by_seq = defaultdict(list)
    for match in matches:
        seq = match['sequence']
        matches_by_seq[seq].append(match)
    unique_matches = []
    for seq, matches in matches_by_seq.items():
        total_area = 0.0
        for match in matches:
            total_area += float(match['attr']['total_area'])
        total_area /= float(len(matches))
        unique_match = peptagram.proteins.new_match(seq)
        unique_match['attr']['total_area'] = total_area
        unique_matches.append(unique_match)
    source['matches'] = unique_matches


def average_across_sources(protein):
    sources = protein['sources']
    area_by_seq = defaultdict(float)
    for source in sources:
        for match in source['matches']:
            seq = match['sequence']
            area_by_seq[seq] += match['attr']['total_area']
    for i_source, source in enumerate(sources):
        for match in source['matches']:
            seq = match['sequence']
            sum_total_area = area_by_seq[seq]
            match['intensity'] = match['attr']['total_area'] / \
                sum_total_area


def convert_skyline_to_peptagram(params, print_fn=sys.stdout.write):
    fname = params['skyline_csv']
    proteins = {}
    for row in csv.DictReader(open(fname)):
        seqid = peptagram.proteins.clean_seqid(row['Protein Name'])
        if seqid not in proteins:
            proteins[seqid] = peptagram.proteins.new_protein(seqid)

    i_source_from_source = {}
    sources = []
    for row in csv.DictReader(open(fname)):
        seqid = peptagram.proteins.clean_seqid(row['Protein Name'])
        seq = row['Peptide Sequence']
        source = row['Replicate Name']
        protein = proteins[seqid]
        if source not in i_source_from_source:
            i = len(i_source_from_source)
            i_source_from_source[source] = i
            if i > 0:
                add_source_to_proteins(proteins)
            sources.append(source)
        i_source = i_source_from_source[source]

        match = peptagram.proteins.new_match(seq)
        match['attr']['total_area'] = row['Total Area']
        protein['sources'][i_source]['matches'].append(match)

    for seqid, protein in proteins.items():
        for source in protein['sources']:
            average_same_peptides_in_source(source)

    for seqid, protein in proteins.items():
        average_across_sources(protein)

    peptagram.proteins.filter_proteins(proteins, params)
  
    peptagram.proteins.make_graphical_comparison_visualisation({
        'title': params['title'],
        'proteins': proteins,
        'source_labels': sources,
        'color_names': ['A=1', 'A=0'],
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

    self.title('Create Skyline QC Peptagram')

    self.push_text("Create Skyline QC Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED PARAMETERS", 20)
    self.push_labeled_param(
        'skyline_csv', 'Skyline .csv', 'skyline.csv', load_file_text='select .csv')
    self.push_labeled_param(
        'fasta', 'Protein sequences in fasta', 'sequences.fasta', load_file_text='select .fasta')
    self.push_labeled_param(
        'title', 'Peptagram title', 'Skyline QC Peptagram')
    self.push_labeled_param(
        'out_dir', 'Output directory', 'peptagram-qc', load_dir_text='select directory')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("OPTIONAL PARAMETERS", 20)

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
    index = convert_skyline_to_peptagram(params, self.print_output)
    callback = lambda: webbrowser.open('file://' + index)
    self.print_output(index, callback)



if __name__ == "__main__":
  if 'test' in sys.argv:
    convert_skyline_to_peptagram(test_params)
  else:
    form = PeptagramForm(800, -50)
    form.mainloop()


