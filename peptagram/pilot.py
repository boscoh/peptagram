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

    for seqid in seqids:
      if seqid not in proteins:
        protein = parse_proteins.new_protein(seqid)
        protein['attr']['seqid'] = seqid
        proteins[seqid] = protein
      else:
        protein = proteins[seqid]

    match_sequence = pilot_match['Sequence']
    match = parse_proteins.new_match(match_sequence)
    match['attr']['missed_cleavages'] = \
        parse_cleavages(pilot_match['Cleavages'])
    match['modifications'] = \
        parse_modifications(
            pilot_match['Modifications'], len(match_sequence))

    def set_attr(pilot_key, attr_key):
      if pilot_key in pilot_match:
        match['attr'][attr_key] = pilot_match[pilot_key]
    set_attr('Sc', 'score')
    set_attr('Score', 'score')
    set_attr('Acq Time', 'retention_time')
    set_attr('Time', 'retention_time')

    protein['sources'][0]['matches'].append(match)

  return proteins



