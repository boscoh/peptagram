 # -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import pymzml


"""
Loads specific MS-MS spectra from an mzML file into 
an existing proteins data structure.
"""

def load_mzml_into_matches(matches, mzml, n_peak=50):
  match_by_scan_id = {}
  for match in matches:
    scan_id = match['attr']['scan_id']
    match_by_scan_id[scan_id] = match

  for spectrum in pymzml.run.Reader(mzml):
    if spectrum['id'] in match_by_scan_id:
      match = match_by_scan_id[spectrum['id']]
      ions = [(mz, i) for mz, i in spectrum.peaks]
      ions.sort(key=lambda i:-i[1])
      match['spectrum'] = ions[:n_peak]


def load_mzml(proteins, i_source, mzml, n_peak=50):
  matches = []
  for protein in proteins.values():
    matches.extend(protein['sources'][i_source]['matches'])
  load_mzml_into_matches(matches, mzml, n_peak)


