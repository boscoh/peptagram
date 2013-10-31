#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) 2013 Bosco Ho
from __future__ import print_function
from pprint import pprint

"""
A peptide mass calculation library.
Given sequences in terms of strings, generates all
the different masses of the peptide depending on the
weighted averages of the different isotypes of each atom.
"""

atom_mass = {
  'H': 1.007825032,
  'O': 15.994914622,
}

ave_atom_mass = {
  'H': 1.00794,
  'O': 15.9994,
}

aa_monoisotopic_mass = {
  'A': 71.037114,
  'B': 110.000000,
  'C': 103.009184,
  'D': 115.026943,
  'E': 129.042593,
  'F': 147.068414,
  'G': 57.021464,
  'H': 137.058912,
  'I': 113.084064,
  'J': 110.000000,
  'K': 128.094963,
  'L': 113.084064,
  'M': 131.040485,
  'N': 114.042927,
  'O': 110.000000,
  'P': 97.052764,
  'Q': 128.058578,
  'R': 156.101111,
  'S': 87.032028,
  'T': 101.047678,
  'U': 110.000000,
  'V': 99.068414,
  'W': 186.079313,
  'X': 110.000000,
  'Y': 163.063329,
  'Z': 110.000000
}

aa_average_mass = {
  'A': 71.0788,
  'B': 110.0000,
  'C': 103.1388,
  'D': 115.0886,
  'E': 129.1155,
  'F': 147.1766,
  'G': 57.0519,
  'H': 137.1411,
  'I': 113.1594,
  'J': 110,
  'K': 128.1741,
  'L': 113.1594,
  'M': 131.1926,
  'N': 114.1038,
  'O': 110,
  'P': 97.1167,
  'Q': 128.1307,
  'R': 156.1875,
  'S': 87.0782,
  'T': 101.1051,
  'U': 110,
  'V': 99.1326,
  'W': 186.2132,
  'X': 110,
  'Y': 163.1760,
  'Z': 110
}

ion_type_props = {
  'b': { 
    'atoms_to_add': [], # Cterm O+ loses electron to form triple bond to C
    'atoms_to_sub': [],
    'is_nterm': True, 
  },
  'y': { 
    'atoms_to_add': ['H', 'H', 'O'], # H+ on Nterm NH, OH on Cterm CO
    'atoms_to_sub': [],
    'is_nterm': False,
  },
}


def parse_aa_masses(sequence, aa_mass, modified_masses=[]):
  masses = [aa_mass[aa] for aa in sequence]
  for modified in modified_masses:
    masses[modified['i']] = modified['mass']
  return masses


def ion_label(i, charge, ion_type):
  if charge == 1:
    return '%s%d' % (ion_type, i)
  return '%s%d(%d+)' % (ion_type, i, charge)


def calculate_neutral_mass(aa_masses, ion_type):
  props = ion_type_props[ion_type]
  mass = sum(aa_masses) 
  mass += sum([atom_mass[a] for a in props['atoms_to_add']])
  mass -= sum([atom_mass[a] for a in props['atoms_to_sub']])
  return mass


def calculate_mz(aa_masses, charge, ion_type):
  charged_mass = calculate_neutral_mass(aa_masses, ion_type)
  charged_mass += atom_mass['H']*charge 
  if charge == 1:
    return charged_mass
  else:
    return charged_mass/charge


def calculate_peaks(aa_masses, charge, ion_type):
  peaks = []
  for i in range(len(aa_masses)):
    if ion_type_props[ion_type]['is_nterm']:
      fragment = aa_masses[:i+1]
    else:
      fragment = aa_masses[i:]
    mz = calculate_mz(fragment, charge, ion_type)
    # round to 4 decimal places
    mz = int(mz*1E4)/1E4
    label = ion_label(len(fragment), charge, ion_type)
    peaks.append((mz, label))
  return peaks


def map_matched_ions(
    ion_type, sequence, peaks, mz_delta=0.8, modified_aa_masses=[],
    aa_mass=aa_monoisotopic_mass):
  aa_masses = parse_aa_masses(sequence, aa_mass, modified_aa_masses)
  ion = ion_type[0]
  pieces = ion_type.split('(')
  charge = pieces[1][0] if len(pieces) > 1 else 1
  theory_peaks = calculate_peaks(aa_masses, charge, ion)
  matched = []
  for theory_mz, label in theory_peaks:
    for mz, intensity in peaks:
      diff = theory_mz - mz
      error = int(diff/theory_mz*1E6)
      if abs(diff) <= mz_delta:
        matched.append([mz, intensity, label, diff, error])
        break
  return matched


if __name__ == "__main__":
  seq = 'MSAFLLTKR'
  spectrum = [[622.494, 2192.0],
   [615.666, 1811.0],
   [619.456, 1810.0],
   [628.245, 1772.0],
   [616.511, 1016.0],
   [614.583, 780.0],
   [629.003, 673.0],
   [490.215, 563.0],
   [606.378, 533.0],
   [563.797, 430.0],
   [610.871, 421.0],
   [472.861, 376.0],
   [607.103, 346.0],
   [584.56, 341.0],
   [605.67, 334.0],
   [620.333, 310.0],
   [375.222, 299.0],
   [625.311, 263.0],
   [488.298, 254.0],
   [786.349, 223.0],
   [899.377, 221.0],
   [506.51, 215.0],
   [881.571, 215.0],
   [601.8, 207.0],
   [514.999, 202.0],
   [598.419, 201.0],
   [608.138, 199.0],
   [593.698, 196.0],
   [759.359, 196.0],
   [530.555, 186.0],
   [623.222, 186.0],
   [602.825, 184.0],
   [542.091, 180.0],
   [261.063, 169.0],
   [627.002, 160.0],
   [913.689, 160.0],
   [592.707, 157.0],
   [825.581, 157.0],
   [471.826, 156.0],
   [505.835, 155.0],
   [609.168, 154.0],
   [768.284, 151.0],
   [563.027, 146.0],
   [570.712, 146.0],
   [580.803, 145.0],
   [979.473, 137.0],
   [596.779, 130.0],
   [541.113, 129.0],
   [777.839, 186.0],
   [550.083, 159.0]];
  delta_mass = 0.5
  modified_aa_masses = [{ 'i':3, 'mass':100.0 }]
  modified_aa_masses = []
  aa_masses = parse_aa_masses(seq, aa_monoisotopic_mass, modified_aa_masses)
  pprint(aa_masses)
  bion_peaks = calculate_peaks(aa_masses, 1, 'b')
  yion_peaks = calculate_peaks(aa_masses, 1, 'y')
  pprint(yion_peaks)
  peaks = []
  matched = map_matched_ions('y', seq, spectrum, delta_mass, modified_aa_masses, aa_monoisotopic_mass)
  matched = map_matched_ions('b', seq, spectrum, delta_mass, modified_aa_masses, aa_monoisotopic_mass)
  pprint(matched);

