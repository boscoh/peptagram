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
        q_good=0.5, 
        q_cutoff=0.75)

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


def make_protein(i_group, protein_group):
    descriptions = protein_group['protein description'].split(' / ')
    coverage_str = str(protein_group['protein sequence coverage (%)'])
    if ';' in coverage_str:
      coverage =  float(get_first(coverage_str, ';'))
    else:
      coverage =  float(get_first(coverage_str, '/'))
    seqs = protein_group['protein sequence'].split('/')
    seqids = [desc.split()[0] for desc in descriptions]
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
    return protein


def make_match(psm, modification_table):
    extracted_peptide_sequence, modifications = parse_peptide(
        psm['peptide sequence'], modification_table)
    peptide_sequence = psm['base peptide sequence']
    if extracted_peptide_sequence != peptide_sequence:
      logger.warning("Peptide sequences don't match: " + psm['peptide sequence'] + " " + extracted_peptide_sequence + " " + peptide_sequence)

    q_value = float(psm['q-value (%)'])

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
      'modifications': [],
      'intensity': 1.0,
      'i': -1,
    }
    if modifications:
      for modification in modifications:
        modification['mass'] = parse.round_decimal(modification['mass'], 4)
      match['modifications'] = modifications
      modified_sequence = psm['peptide sequence'].split('.')[1]
      match['attr']['modified_sequence'] = modified_sequence

    return match


def get_i_source(proteins, sources, source):
  for i_test_source, test_source in enumerate(sources):
    if source == test_source:
      i_source = i_test_source
      break
  else:
    i_source = len(sources)
    sources.append(source)
    if i_source > 0:
      for protein in proteins.values():
        protein['sources'].append({ 'matches':[] })
  return i_source


def get_proteins_and_sources(
      protein_groups_fname, 
      psm_fname, 
      modifications_fname=None,
      q_good=0.0, 
      q_cutoff=0.75):

  is_debug = logger.root.level <= logging.DEBUG

  dump_dir = os.path.dirname(protein_groups_fname)

  modification_table = {}
  if modifications_fname is not None:
    modification_table = read_modification_dict(modifications_fname)
  
  proteins = {}
  dict_dump_writer = parse.DictListWriter(is_debug, os.path.join(dump_dir, 'protein_groups.dump'))
  for i_group, protein_group in enumerate(parse.read_tsv(protein_groups_fname)):
    protein = make_protein(i_group, protein_group)
    proteins[protein['attr']['seqid']] = protein
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
  i_source_from_source = {}
  sources = []
  for psm in parse.read_tsv(psm_fname):
    dict_dump_writer.dump_dict(psm)

    match = make_match(psm, modification_table)
    match['intensity']  = parse_proteins.calc_intensity(
       match['attr']['q_value'], q_good, q_cutoff)
    if match['attr']['q_value'] > q_cutoff:
      continue
    peptide_sequence = match['sequence']

    n_match += 1

    protein = None
    descriptions = psm['protein description'].split(' / ')
    peptide_seqids = [d.split()[0] for d in descriptions]
    for peptide_seqid in peptide_seqids:
      if peptide_seqid in protein_by_seqid:
        test_protein = protein_by_seqid[peptide_seqid]
        sequence = protein_by_seqid[peptide_seqid]['sequence']
        if peptide_sequence in sequence:
          protein = test_protein
          break
    else:
      logger.warning("Couldn't find protein for %s" % (peptide_sequence))
      continue
    match['i'] = sequence.find(peptide_sequence)

    n_match_assigned += 1
    i_source = get_i_source(proteins, sources, psm['filename'])
    protein['sources'][i_source]['matches'].append(match)

  dict_dump_writer.close()

  dump = os.path.join(dump_dir, 'proteins.dump')
  if logger.root.level <= logging.DEBUG:
    logger.debug('Dumping proteins data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  logger.info("Assigned {}/{} of PSMs.tsv to protein_groups.tsv".format(n_match_assigned, n_match))

  return proteins, sources


