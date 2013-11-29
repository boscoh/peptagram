 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import pymzml

def load_mzml(proteins, i_source, mzml, n_peak=50):
  peptide_by_scan_id = {}
  for protein in proteins.values():
    source = protein['sources'][i_source]
    for peptide in source['peptides']:
      scan_id = peptide['attr']['scan_id']
      peptide_by_scan_id[scan_id] = peptide
  for spectrum in pymzml.run.Reader(mzml):
    if spectrum['id'] in peptide_by_scan_id:
      peptide = peptide_by_scan_id[spectrum['id']]
      ions = [(mz, i) for mz, i in spectrum.peaks]
      ions.sort(key=lambda i:-i[1])
      peptide['spectrum'] = ions[:n_peak]


