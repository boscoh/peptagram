 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import pyteomics.mzml
import parse


def make_peaks(mz_str, intensity_str, n_peak):
  x_vals = map(float, mz_str)
  y_vals = map(float, intensity_str)
  ions = [(x, y) for x, y in zip(x_vals, y_vals)]
  ions.sort(key=lambda i:-i[1])
  ions = ions[:n_peak]
  return ions


def load_mzml(proteins, i_source, mzml, n_peak=50):
  peptide_by_scan_id = {}
  for protein in proteins.values():
    source = protein['sources'][i_source]
    for peptide in source['peptides']:
      scan_id = peptide['attr']['scan_id']
      peptide_by_scan_id[scan_id] = peptide
  for spectrum in pyteomics.mzml.read(mzml):
    scan_id = int(spectrum['id'].split()[-1].split('=')[-1])
    peptide['attr']['scan_id'] = scan_id
    if scan_id in peptide_by_scan_id:
      peptide = peptide_by_scan_id[scan_id]
      peaks = make_peaks(
          spectrum['m/z array'], spectrum['intensity array'], n_peak)
      peptide['spectrum'] = peaks
      # if 'precursorList' in spectrum:
      #   mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
      #   mz = float('%.5f' % mz)
      #   peptide['attr']['m/z mzml'] = parse.round_decimal(mz, 4)


