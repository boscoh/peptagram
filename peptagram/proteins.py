 # -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import math
import os
import json
import copy
import glob
import shutil
from pprint import pprint

import logging
logger = logging.getLogger('proteins')

import uniprot
import parse


"""
Data structure manipulation routines for the key protein dictionary
created with all the different parsers. 

The protein dictionary is used to collect peptide-spectrum-matches
organized around protein sequences. Routines have been provided
to filter and manipulate this dictionary.

The protein dictionary is then processed to generate peptatgram
visualizations.
"""


this_dir = os.path.abspath(os.path.dirname(__file__))


def clean_seqid(seqid):
  if '|' in seqid:
    return seqid.split('|')[1]
  else:
    return seqid


def new_protein(seqid):
  return {
    'attr': { 'seqid':seqid, 'other_seqids':[] },
    'description': '',
    'sequence': '',
    'sources': [{'matches': [], }],
  }


def new_match(peptide_sequence):
  return {
    'sequence': peptide_sequence,
    'intensity': 1.0,
    'modifications': [],
    'mask': 0.0,
    'attr': {}
  }


def match_iterator(proteins):
  for seqid in proteins:
    protein = proteins[seqid]
    for source in protein['sources']:
      for match in source['matches']:
        yield seqid, match


def do_matches(proteins, do_match_fn):
  for protein in proteins.values():
    for source in protein['sources']:
      for match in source['matches']:
        do_match_fn(match)


def determine_unique_matches(proteins):
  seqids_by_seq = {}
  for seqid, match in match_iterator(proteins):
    seq = match['sequence']
    if seq not in seqids_by_seq:
      seqids_by_seq[seq] = set()
    seqids_by_seq[seq].add(seqid)
  for seqid, match in match_iterator(proteins):
    seq = match['sequence']
    n_seqid = len(seqids_by_seq[seq])
    match['attr']['is_unique'] = (n_seqid == 1)


def check_missing_fields(proteins):
  for seqid, match in match_iterator(proteins):
    if 'modifications' not in match:
      match['modifications'] = []
    if 'intensity' not in match:
      match['intensity'] = 1
    if 'mask' not in match:
      match['mask'] = 0.0


def count_matches(proteins):
  for seqid in proteins.keys():
    protein = proteins[seqid]
    n_match = 0
    n_slice_populated = 0
    n_unique_match = 0
    for source in protein['sources']:
      matches = source['matches']
      unique_matches = [p for p in matches if p['attr']['is_unique']]
      n_unique_match += len(unique_matches)
      n_match += len(matches)
      if len(matches) > 0:
        n_slice_populated += 1
    protein['attr']['n_slice_populated'] = n_slice_populated
    protein['attr']['n_match_unique'] = n_unique_match
    protein['attr']['n_match'] = n_match


def delete_empty_proteins(proteins):
  for seqid in proteins.keys():
    protein = proteins[seqid]
    n_match = 0
    for source in protein['sources']:
      n_match += len(source['matches'])
    if n_match == 0:
      del proteins[seqid]


def delete_spectrum(match):
  if 'spectrum' in match:
    del match['spectrum']


def delete_matches(proteins, is_deleteable_fn):
  for protein in proteins.values():
    for source in protein['sources']:
      matches = source['matches']
      n_match = len(matches)
      for i_match in reversed(range(n_match)):
        if is_deleteable_fn(matches[i_match]):
          del matches[i_match]
  delete_empty_proteins(proteins)


def is_missed_cleavage(match):
  if 'missed_cleavages' in match['attr'] and \
      match['attr']['missed_cleavages'] > 0:
    return True
  return False


def is_tryptic(match):
  if 'missed_cleavages' not in match['attr']:
    return True
  return match['attr']['missed_cleavages'] == 0


def is_modified_peptide(match):
  if 'modifications' not in match:
    return True
  if len(match['modifications']) == 0:
    return True
  return False


def calc_intensity(x, high, low):
  return (x-low)/(high-low)*.8 + .2


def calc_minus_log_intensity(x, high, low):
  return calc_intensity(
    -math.log(x), 
    -math.log(high), 
    -math.log(low))


def find_peptide_positions_in_proteins(proteins):
  protein_list = proteins.values()
  if len(protein_list) == 0:
    return
  n_source = len(protein_list[0]['sources'])
  seqids_by_sequence = {}
  for seqid in proteins:
    protein = proteins[seqid]
    for i_source in range(n_source):
      source = protein['sources'][i_source]
      matches = source['matches']
      for match in matches:
        sequence = match['sequence']
        if sequence not in seqids_by_sequence:
          seqids_by_sequence[sequence] = set()
        seqids_by_sequence[sequence].add(seqid)
  for seqid in proteins:
    protein = proteins[seqid]
    for i_source in range(n_source):
      protein = proteins[seqid]
      source = protein['sources'][i_source]
      matches = source['matches']
      for match in matches:
        sequence = match['sequence']
        for test_seqid in seqids_by_sequence[sequence]:
          if seqid != test_seqid:
            if 'other_seqids' not in match['attr']:
              match['attr']['other_seqids'] = []
            match['attr']['other_seqids'].append(test_seqid)


def change_seqids_in_proteins(proteins, clean_seqid):
  seqids = proteins.keys()
  for seqid in seqids:
    new_seqid = clean_seqid(seqid)
    if 'attr' in proteins[seqid]:
      if 'other_seqids' in proteins[seqid]['attr']:
        other_seqids = proteins[seqid]['attr']['other_seqids']
        other_seqids = map(clean_seqid, other_seqids)
        proteins[seqid]['attr']['other_seqids'] = other_seqids
    if seqid != new_seqid:
      proteins[new_seqid] = proteins[seqid]
      del proteins[seqid]
      if 'attr' in proteins[new_seqid]:
        proteins[new_seqid]['attr']['seqid'] = new_seqid


def calculate_peptide_positions(proteins, iso_leu_isomerism=False):
  for seqid in proteins:
    protein = proteins[seqid]
    protein_sequence = protein['sequence']
    if iso_leu_isomerism:
      protein_sequence = protein_sequence.replace("L", "I")
    for source in protein['sources']:
      matches = source['matches']
      for i_match in reversed(range(len(matches))):
        match = matches[i_match]
        peptide_sequence = match['sequence']
        if iso_leu_isomerism:
          peptide_sequence = peptide_sequence.replace("L", "I")
        i = protein_sequence.find(peptide_sequence)
        if i < 0:
          logger.debug("'{}' not found in {}".format(peptide_sequence, seqid))
          del matches[i_match]
          continue
        match['i'] = i 
      for match in matches:
        if 'i' not in match:
          pprint(match)
          raise ValueError


def load_fastas_into_proteins(
    proteins, fastas, clean_seqid=None, iso_leu_isomerism=False):
  if clean_seqid:
    change_seqids_in_proteins(proteins, clean_seqid)
    change_seqids_in_proteins(fastas, clean_seqid)
  for seqid in proteins.keys():
    protein = proteins[seqid]
    if seqid not in fastas:
      logger.debug("%s of proteins not found in fastas" % seqid)
      del proteins[seqid]
      continue
    protein_sequence = fastas[seqid]['sequence']
    protein['description'] = fastas[seqid]['description']
    protein['sequence'] = protein_sequence
    protein['attr']['length'] = len(protein_sequence)
  calculate_peptide_positions(proteins, iso_leu_isomerism)


def load_sequences_from_uniprot(proteins, clean_seqid=None, cache_basename=None):
  if clean_seqid:
    change_seqids_in_proteins(proteins, clean_seqid)
  seqids = []
  for seqid in proteins:
    seqids.append(seqid)
    if 'other_seqids' in proteins[seqid]['attr']:
      seqids.extend(proteins[seqid]['attr']['other_seqids'])
  uniprot_data = uniprot.get_metadata_with_some_seqid_conversions(seqids, cache_basename)
  load_fastas_into_proteins(proteins, uniprot_data)
  if cache_basename:
    uniprot.write_fasta(cache_basename+'.fasta', uniprot_data, uniprot_data.keys())


def load_fasta_db_into_proteins(
    proteins, fasta_db, clean_seqid=None, iso_leu_isomerism=False):
  seqids, fastas = parse.read_fasta(fasta_db)
  load_fastas_into_proteins(proteins, fastas, clean_seqid, iso_leu_isomerism)


def merge_two_proteins(proteins1, proteins2):
  """
  Merges two proteins structures. In particular, it grafts the
  'sources' together, treating the sources in each proteins as 
  distinct, and maintaing the order.
  """
  if len(proteins1) == 0:
    return proteins2
  seqid = proteins1.keys()[0]
  n_source1 = len(proteins1[seqid]['sources'])
  seqid = proteins2.keys()[0]
  n_source2 = len(proteins2[seqid]['sources'])
  for seqid in proteins2:
    protein2 = proteins2[seqid]
    sources2 = protein2['sources']
    if seqid in proteins1:
      protein1 = proteins1[seqid]
      protein1['sources'].extend(sources2)
    else:
      sources1 = [{'matches': []} for i in range(n_source1)]
      combined_sources = sources1 + sources2
      proteins1[seqid] = protein2
      proteins1[seqid]['sources'] = combined_sources
  for seqid in proteins1:
    matches = proteins1[seqid]['sources']
    if len(matches) == n_source1:
      matches.extend([{'matches': []} for i in range(n_source2)])
  return proteins1


def save_data_js(data, js_fname):
  f = open(js_fname, 'w')
  f.write('var data = \n')
  f.write(json.dumps(data, indent=None))
  f.close()


def transfer_files(in_dir, out_dir):
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
  for src in glob.glob(os.path.join(in_dir, '*')):
    dst = os.path.join(out_dir, os.path.basename(src))
    shutil.copy(src, dst)



def do_peptides(proteins, peptide_fn):
  for protein in proteins.values():
    for source in protein['sources']:
      for match in source['matches']:
        peptide_fn(match)



aa2res = {  '<': 'NME', '>': 'ACE', 'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E':
'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T':
'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}



def mod_str(peptide):
  s = ''
  if 'modifications' not in peptide:
    return 
  for mod in peptide['modifications']:
    i_peptide = mod['i']
    if i_peptide < 0:
      res = 'NTR'
    elif i_peptide >= len(peptide['sequence']):
      res = 'CTR'
    else:
      aa = peptide['sequence'][i_peptide]
      if aa in aa2res:
        res = aa2res[aa]
      else:
        res = 'XXX'
    i_protein = peptide['i'] + i_peptide
    s += "<br>&nbsp; %s-%d(%d)" % (
        res,
        i_protein + 1,
        i_peptide + 1)
    if 'mass' in mod:
      s += " M=%.2f" % float(mod['mass'])
    if 'type' in mod:
      s += " %s" % mod['type']
  peptide['attr']['modifications'] = s


def filter_proteins(proteins, params):
  """
  Filters protein entries, for use with the GUI. Namely
  params is a dictionary describing potential filtering.
  """
  if 'include_seqids' in params and params['include_seqids']:
    seqids = parse.read_word_file(params['include_seqids'])
    for seqid in proteins:
      if seqid not in seqids:
        del proteins[seqid]

  if 'exculde_seqids' in params and params['exclude_seqids']:
    for seqid in parse.read_word_file(params['exclude_seqids']):
      if seqid in proteins:
        del proteins[seqid]

  if 'fasta' in params and params['fasta']:
    load_fasta_db_into_proteins(
        proteins, params['fasta'], clean_seqid=clean_seqid)

  if 'include_msms' in params and params['include_msms'] == 0:
    do_matches(proteins, delete_spectrum)

  if 'match_filter'in params:
    if params['match_filter'] == 1: # tryptic
      delete_matches(proteins, is_missed_cleavage)
    elif params['match_filter'] == 2: # semitryptic
      delete_matches(proteins, is_tryptic)
    elif params['match_filter'] == 3: # modified
      delete_matches(proteins, is_modified_peptide)


def make_graphical_comparison_visualisation(data, out_dir=None):
  # sanity checks
  proteins = data['proteins']
  determine_unique_matches(proteins)
  delete_empty_proteins(proteins)
  check_missing_fields(proteins)
  count_matches(proteins)
  do_peptides(proteins, mod_str)

  find_peptide_positions_in_proteins(proteins)
  for seqid, protein in proteins.items():
    for source in protein['sources']:
      matches = source['matches']
      matches.sort(key=lambda match: len(match['sequence']))
      matches.sort(key=lambda match: match['i'])

  if 'source_labels' not in data:
    data['source_labels'] = []
  if 'color_names' not in data:
    data['color_names'] = ['', '', '']
  if 'mask_labels' not in data:
    data['mask_labels'] = []

  if out_dir is None and 'out_dir' in data:
    out_dir = data['out_dir']
    data = data.copy()
    del data['out_dir']
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  save_data_js(data, os.path.join(out_dir, 'data.js'))
  transfer_files(os.path.join(this_dir, 'templates/comparison'), out_dir)
  transfer_files(os.path.join(this_dir, 'templates/js'), os.path.join(out_dir, 'js'))
  index_html = os.path.abspath(os.path.join(out_dir, 'index.html'))
  logger.info('Made peptograph in "' + index_html + '"')


def make_sequence_overview_visualisation(data, out_dir=None):
  # sanity checks
  proteins = data['proteins']
  determine_unique_matches(proteins)
  delete_empty_proteins(proteins)
  check_missing_fields(proteins)
  find_peptide_positions_in_proteins(proteins)
  for seqid, protein in proteins.items():
    for source in protein['sources']:
      matches = source['matches']
      matches.sort(key=lambda match: len(match['sequence']))
      matches.sort(key=lambda match: match['i'])
  if 'source_labels' not in data:
    data['source_labels'] = []
  if 'color_names' not in data:
    data['color_names'] = ['', '', '']
  if 'mask_labels' not in data:
    data['mask_labels'] = []

  if out_dir is None and 'out_dir' in data:
    out_dir = data['out_dir']
    data = data.copy()
    del data['out_dir']
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  save_data_js(data, os.path.join(out_dir, 'data.js'))
  transfer_files(os.path.join(this_dir, 'templates/overview'), out_dir)
  transfer_files(os.path.join(this_dir, 'templates/js'), os.path.join(out_dir, 'js'))
  index_html = os.path.abspath(os.path.join(out_dir, 'index.html'))
  logger.info('Made peptograph in "' + index_html + '"')





