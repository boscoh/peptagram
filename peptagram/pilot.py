 # -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import csv

import logging
logger = logging.getLogger('pilot')

import proteins as parse_proteins


"""
Parser for AB/X Protein Pilot search results saved as CSV text file

Main API entry:

  get_proteins(fname)

  returns a dictionary that organizes peptide-spectrum-matches
  around proteins.

"""


def get_delimiter(fname):
  if fname.endswith('csv'):
    return ','
  else:
    return '\t'


def read_matches(pilot_summary):
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
    mod_type, position = piece.split('@')
    if position == 'N-term': 
      i = -1
    elif position == 'C-term':
      i = len_seq
    else:
      i = int(position) - 1
    modification = { 'i':i, 'type':mod_type }
    modifications.append(modification)
  return modifications
  

def get_proteins(fname):
  proteins = {}

  for pilot_match in read_matches(fname):

    seqids = map(
      parse_proteins.clean_seqid,
      pilot_match['Accessions'].split(';'))

    seqid = seqids[0]
    if seqid not in proteins:
      protein = parse_proteins.new_protein(seqid)
      proteins[seqid] = protein
      protein['attr']['other_seqids'] = seqids[1:]
      protein['attr']['seqid'] = seqid
    else:
      protein = proteins[seqid]

    match_sequence = pilot_match['Sequence']
    match = parse_proteins.new_match(match_sequence)
    match['attr']['missed_cleavages'] = \
        parse_cleavages(pilot_match['Cleavages'])
    match['modifications'] = \
        parse_modifications(
            pilot_match['Modifications'], len(match_sequence))

    if 'Sc' in pilot_match:
      match['attr']['score'] = pilot_match['Sc']
    elif 'Score' in pilot_match:
      match['attr']['score'] = pilot_match['Score']
    if 'Acq Time' in pilot_match:
      match['attr']['retention_time'] = pilot_match['Acq Time']
    elif 'Time' in pilot_match:
      match['attr']['retention_time'] = pilot_match['Time']
    match['attr']

    protein['sources'][0]['matches'].append(match)

  return proteins


if __name__ == "__main__":
  logging.basicConfig(level=logging.INFO)

