 # -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import os
import json
import csv

import logging
logger = logging.getLogger('morpheus')

import parse
import proteins as parse_proteins
import peptidemass


"""
Parser for Morpheus .TSV files for mass-spec search results

Main API entry:

  get_proteins(
        protein_groups_fname, 
        psm_fname, 
        modifications_fname=None,
        q_okay=0.5, q_cutoff=0.75)

  returns a dictionary that organizes peptide-spectrum-matches
  around proteins.
"""


def read_modification_dict(modifications_tsv):
  result = {}
  for entry in parse.read_tsv(modifications_tsv):
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
          logger.debug('Warning: modification {} unknown'.format(mod_str))
          continue
        modification = {
            'i':i, 
            'mass': float(modification_dict[mod_str]),
            'type': mod_str,
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


def get_proteins(
      protein_groups_fname, 
      psm_fname, 
      modifications_fname=None,
      q_okay=0.5, q_cutoff=0.75):

  is_debug = logger.root.level <= logging.DEBUG

  dump_dir = os.path.dirname(protein_groups_fname)

  if modifications_fname is not None:
    modification_table = read_modification_dict(modifications_fname)
  else:
    modification_table = {}
  
  proteins = {}
  dict_dump_writer = parse.DictListWriter(is_debug, os.path.join(dump_dir, 'protein_groups.dump'))
  for i_group, protein_group in enumerate(parse.read_tsv(protein_groups_fname)):
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
      'sources': [{ 'matches':[] }]
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

  dict_dump_writer = parse.DictListWriter(is_debug, os.path.join(dump_dir, 'peptides.dump'))
  n_match = 0
  n_match_assigned = 0

  for psm in parse.read_tsv(psm_fname):

    dict_dump_writer.dump_dict(psm)
    descriptions = psm['protein description'].split(' / ')
    extracted_peptide_sequence, modifications = parse_peptide(
        psm['peptide sequence'], modification_table)
    peptide_sequence = psm['base peptide sequence']
    if extracted_peptide_sequence != peptide_sequence:
      logger.warning("Peptide sequences don't match: " + psm['peptide sequence'] + " " + extracted_peptide_sequence + " " + peptide_sequence)

    protein = None
    peptide_seqids = [d.split()[0] for d in descriptions]
    for peptide_seqid in peptide_seqids:
      if peptide_seqid in protein_by_seqid:
        test_protein = protein_by_seqid[peptide_seqid]
        sequence = test_protein['sequence']
        i = sequence.find(peptide_sequence)
        if i < 0:
          continue
        else:
          protein = test_protein
        break
    if protein is None:
      continue

    n_match += 1
    n_match_assigned += 1

    q_value = float(psm['q-value (%)'])
    if q_value > q_cutoff:
      continue

    intensity = parse_proteins.calc_intensity(
       q_value, 0.0, q_cutoff)

    mask = 0
    if q_value < q_okay:
      mask = 1

    if 'scan number' in psm:
      scan_id = psm['scan number']
    elif 'spectrum number' in psm:
      scan_id = psm['spectrum number']
    else:
      scan_id = ''
    if 'retention time (min)' in psm:
      time = parse.round_decimal(psm['retention time (min)'], 4)
    elif 'retention time (minutes)' in psm:
      time = parse.round_decimal(psm['retention time (minutes)'], 4)
    else:
      time = ''

    match = {
      'sequence': peptide_sequence,
      'attr': {
        'scan_id': scan_id, 
        'retention_time': time,
        'morpheus_score': parse.round_decimal(psm['morpheus score'], 4),
        'mass': parse.round_decimal(psm['precursor mass (da)'], 4),
        'mass_diff': parse.round_decimal(psm['precursor mass error (da)'], 4),
        'm/z': parse.round_decimal(psm['precursor m/z'], 4),
        'source': parse.basename(psm['filename']),
        'missed_cleavages': int(psm['missed cleavages']),
        'q_value': q_value,
      },
      'intensity': intensity,
      'mask': mask,
      'i': i,
    }
    if modifications:
      for modification in modifications:
        modification['mass'] = parse.round_decimal(modification['mass'], 4)
      match['modifications'] = modifications
      modified_sequence = psm['peptide sequence'].split('.')[1]
      match['attr']['modified_sequence'] = modified_sequence

    protein['sources'][0]['matches'].append(match)

  dict_dump_writer.close()

  dump = os.path.join(dump_dir, 'proteins.dump')
  if logger.root.level <= logging.DEBUG:
    logger.debug('Dumping proteins data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  logger.info("Assigned {}/{} of PSMs.tsv to protein_groups.tsv".format(n_match_assigned, n_match))

  return proteins


if __name__ == '__main__':
  logging.basicConfig(level=logging.DEBUG)
  proteins = get_proteins(
      '../example/morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
      '../example/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
      '../example/morpheus/modifications.tsv')


