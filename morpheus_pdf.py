 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import logging
logging.basicConfig(level=logging.WARNING)
import os
import datetime

import peptagram.morpheus
import peptagram.mzml
import peptagram.proteins
import peptagram.pdf


protein_group = 'output/2010-12-05_Yeast_Trypsin_FT_HCD_Rep3.protein_groups.tsv'
psm = 'output/2010-12-05_Yeast_Trypsin_FT_HCD_Rep3.PSMs.tsv'
modificiations = 'modifications.tsv'
mzml = '2010-12-05_Yeast_Trypsin_FT_HCD_Rep3.mzML'
pdf = 'modified_peptides.pdf'
title = "Morpheus Search Results for Modified Peptides: 2010-12-05_Yeast_Trypsin_FT_HCD_Rep3"
author = "Craig Wenger"
q_okay = 0.1 # for coloring purposes in peptagram
q_cutoff = 0.5
n_peak = 400
mz_delta = 0.01
metadata_keys = [
    'modified_sequence', 'modifications', 'source', 
    'scan_id', 'retention_time', 'm/z', 'mass', 
    'mass_diff', 'morpheus_score', 'q_value']

print("Reading Morpheus...")
proteins = peptagram.morpheus.get_proteins(
    protein_group, psm, modificiations, 
    q_okay=q_okay, q_cutoff=q_cutoff)

print("Reading mzml...")
i_source = 0
peptagram.mzml.load_mzml(
    proteins, i_source, mzml, n_peak)

sorted_matches = []
for seqid, protein in proteins.items():
  for source in protein['sources']:
    for match in source['matches']:
      if 'modifications' in match and match['modifications']:
        scan_id = match['attr']['scan_id']
        sorted_matches.append((scan_id, match))
sorted_matches.sort()

print("Making pdf...")
doc = peptagram.pdf.PeptagramDoc(
    pdf, title, author, n_peak, mz_delta, metadata_keys)
for scan_id, match in sorted_matches:
    doc.add_match_page(match)
doc.build()

os.system('open ' + pdf)



