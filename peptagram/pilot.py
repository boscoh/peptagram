 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import csv
import logging

logger = logging.getLogger('pilot')

import proteins as parse_proteins


"""
Parser for AB/X peptide results from Protein Pilot
"""


def get_delimiter(fname):
  if fname.endswith('csv'):
    return ','
  else:
    return '\t'


def read_peptides(pilot_summary):
  delimiter = get_delimiter(pilot_summary)
  with open(pilot_summary, 'Ur') as f:
    for entry in csv.DictReader(f, delimiter=delimiter):
      yield entry


def read_headers(pilot_summary):
  delimiter = get_delimiter(pilot_summary)
  with open(pilot_summary, 'Ur') as f:
    reader = csv.reader(f, delimiter=delimiter)
    field_names = reader.next()
  return field_names


# LGSVVNNEENTCSDKRMKPFEEGHGITQVDK: cleavages: missed K-R@15; missed R-M@16
def parse_cleavages(s):
  return len(s.split(';'))


# LGSVVNNEENTCSDK: modifications: Carbamidomethyl(C)@12; Sulfide(D)@14
def parse_modifications(s, len_seq):
  if not s:
    return []
  modifications = []
  for piece in s.split(';'):
    l, r = piece.split('@')
    if r == 'N-term': 
      i = -1
    elif r == 'C-term':
      i = len_seq
    else:
      i = int(r) - 1
    modification = { 'i':i }
    modifications.append(modification)
  return modifications
  

def get_proteins(fname):
  proteins = {}
  print(read_headers(fname))
  for peptide in read_peptides(fname):
    seqid = parse_proteins.clean_seqid(peptide['Accessions'])
    if seqid not in proteins:
      protein = parse_proteins.new_protein(seqid)
      proteins[seqid] = protein
    else:
      protein = proteins[seqid]
    peptide_sequence = peptide['Sequence']
    match = parse_proteins.new_match(peptide_sequence)
    match['attr']['missed_cleavages'] = \
        parse_cleavages(peptide['Cleavages'])
    match['modifications'] = \
        parse_modifications(peptide['Modifications'], len(peptide_sequence))
    protein['sources'][0]['matches'].append(match)
  return proteins


if __name__ == "__main__":
  logging.basicConfig(level=logging.INFO)

