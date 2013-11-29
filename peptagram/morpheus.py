 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import os
import json
import csv
import logging
import math

import parse
import proteins as parse_proteins
import peptidemass

"""
Parser for Morpheus files
"""


logger = logging.getLogger('morpheus')


def read_tsv_iter(tsv_txt):
  """
  Reads a title top TSV file, converts all keys to lower case
  and returns a list of dictionaries.
  """
  f = open(tsv_txt, "UR")
  titles = [w.lower() for w in parse.split_tab(f.readline())]
  results = []
  for line in f:
    group = {}
    for key, val in zip(titles, parse.split_tab(line)):
      group[key] = parse.parse_string(val)
    yield group
    del group


def read_modification_dict(modifications_tsv):
  result = {}
  for entry in parse.read_tsv(modifications_tsv):
    # pprint(entry)
    key = entry['description']
    mass = float(entry['monoisotopic mass shift (da)'])
    if 'residue' in entry:
      aa = entry['residue']
    elif 'amino acid residue' in entry:
      aa = entry['amino acid residue']
    else:
      aa = 'n/a'
    if aa != 'n/a':
      mass += peptidemass.aa_monoisotopic_mass[aa]
    result[key] = mass
  return result


def parse_peptide(text, modification_dict):
  seq = text.split('.')[1]
  chars = []
  modifications = []
  in_modification = False
  mod_str = ''
  for c in seq:
    if in_modification:
      if c == ']' or c == ')':
        in_modification = False
        i = len(chars) - 1
        if mod_str not in modification_dict:
          log.debug('Warning: modification {} unknown'.format(mod_str))
          continue
        modification = {
            'i':i, 
            'mass': modification_dict[mod_str],
        }
        modifications.append(modification)
      else:
        mod_str += c
    else:
      if c == '[' or c == '(':
        in_modification = True
        mod_str = ''
      else:
        chars.append(c)
  return ''.join(chars), modifications


def get_first(s, delimiter='/'):
  if not isinstance(s, str):
    return s
  elif delimiter not in s:
    return s
  return s.split(delimiter)[0].strip()


class DictListWriter():
  def __init__(self, is_debug, fname):
    self.fname = fname
    self.is_debug = is_debug
    if self.is_debug:
      logger.debug('Dumping dict list to ' + self.fname)
      self.file = open(self.fname, 'w')
      self.file.write('[\n')

  def dump_dict(self, data_dict):
    if self.is_debug:
      pprint(data_dict, stream=self.file)
      self.file.write(',\n')

  def close(self):
    if self.is_debug:
      self.file.write(']\n')
      self.file.close()


def get_proteins(protein_groups_fname, psm_fname, modifications_fname=None):
  is_debug = logger.root.level <= logging.DEBUG
  dump_dir = os.path.dirname(protein_groups_fname)

  if modifications_fname is not None:
    modification_table = read_modification_dict(modifications_fname)
  else:
    modification_table = {}
  
  proteins = {}
  dict_dump_writer = DictListWriter(is_debug, os.path.join(dump_dir, 'protein_groups.dump'))
  for i_group, protein_group in enumerate(read_tsv_iter(protein_groups_fname)):
    descriptions = protein_group['protein description'].split(' / ')
    coverage_str = str(protein_group['protein sequence coverage (%)'])
    if ';' in coverage_str:
      coverage =  float(get_first(coverage_str, ';'))
    else:
      coverage =  float(get_first(coverage_str, '/'))
    seqs = protein_group['protein sequence'].split('/')
    seqids = [desc.split()[0] for desc in descriptions]
    for seqid in seqids:
      if seqid in proteins:
        logger.warning("Different protein groups claim same first seqid", seqid)
    protein = {
      'description': descriptions[0],
      'sequence': seqs[0],
      'other_sequences': seqs[1:],
      'attr': {
        'coverage': parse.round_decimal(coverage, 4),
        'morpheus-score': parse.round_decimal(protein_group['summed morpheus score'], 4),
        'i_group': i_group,
        'other_seqids': seqids[1:],
        'seqid': seqids[0],
      },
      'sources': [{ 'peptides':[] }]
    }
    proteins[seqids[0]] = protein
    dict_dump_writer.dump_dict(protein_group)
  dict_dump_writer.close()

  protein_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    protein_by_seqid[seqid] = protein
    for alt_seqid in protein['attr']['other_seqids']:
      protein_by_seqid[alt_seqid] = protein

  dict_dump_writer = DictListWriter(is_debug, os.path.join(dump_dir, 'peptides.dump'))
  n_peptide = 0
  n_peptide_matched = 0
  for src_peptide in read_tsv_iter(psm_fname):
    dict_dump_writer.dump_dict(src_peptide)
    descriptions = src_peptide['protein description'].split(' / ')
    peptide_seqids = [d.split()[0] for d in descriptions]
    protein = None
    for peptide_seqid in peptide_seqids:
      if peptide_seqid in protein_by_seqid:
        protein = protein_by_seqid[peptide_seqid]
        break
    n_peptide += 1
    if protein is None:
      continue
    n_peptide_matched += 1
    sequence = protein['sequence']
    extracted_peptide_sequence, modifications = parse_peptide(
        src_peptide['peptide sequence'],
        modification_table)
    peptide_sequence = src_peptide['base peptide sequence']
    if extracted_peptide_sequence != peptide_sequence:
      logger.warning("Peptide sequences don't match: " + src_peptide['peptide sequence'] + " " + extracted_peptide_sequence + " " + peptide_sequence)
    i = sequence.find(peptide_sequence)
    if i < 0:
      logger.warning(peptide_sequence + ' not found in ' + protein['attr']['seqid'])
      continue
    q_value = float(src_peptide['q-value (%)'])
    if 'scan number' in src_peptide:
      scan_id = src_peptide['scan number']
    elif 'spectrum number' in src_peptide:
      scan_id = src_peptide['spectrum number']
    else:
      scan_id = ''
    if 'retention time (min)' in src_peptide:
      time = parse.round_decimal(src_peptide['retention time (min)'], 4)
    elif 'retention time (minutes)' in src_peptide:
      time = parse.round_decimal(src_peptide['retention time (minutes)'], 4)
    else:
      time = ''

    peptide = {
      'sequence': peptide_sequence,
      'attr': {
        'scan_id': scan_id, 
        'retention_time': time,
        'morpheus_score': parse.round_decimal(src_peptide['morpheus score'], 4),
        'mass': parse.round_decimal(src_peptide['precursor mass (da)'], 4),
        'mass_diff': parse.round_decimal(src_peptide['precursor mass error (da)'], 4),
        'm/z': parse.round_decimal(src_peptide['precursor m/z'], 4),
        'source': parse.basename(src_peptide['filename']),
        'q_value': q_value,
      },
      'intensity': 1.0 - q_value/100.0,
      'i': i,
    }
    if modifications:
      for modification in modifications:
        modification['mass'] = parse.round_decimal(modification['mass'], 4)
      peptide['attr']['modifications'] = modifications

    protein['sources'][0]['peptides'].append(peptide)

  dict_dump_writer.close()

  dump = os.path.join(dump_dir, 'proteins.dump')
  if logger.root.level <= logging.DEBUG:
    logger.debug('Dumping proteins data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  logger.info("Assigned {}/{} of PSMs.tsv to protein_groups.tsv".format(n_peptide_matched, n_peptide))

  return proteins


if __name__ == '__main__':
  logging.basicConfig(level=logging.DEBUG)
  proteins = get_proteins(
      '../example/morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
      '../example/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
      '../example/morpheus/modifications.tsv')


