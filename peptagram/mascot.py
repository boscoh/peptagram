# -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import math
import json
import urllib
import logging


from parse import parse_string, save_data_dict


logger = logging.getLogger('mascot')


"""
Parser for Mascot dat files 

Some useful formulae:
 - identity_score = -10*Math.log10(1/qmatch)
 - qplughole = homology_score
 - significant if identity_score > homology_score

"""


def piece(s, splitter, i, j=None):
  if j is None:
    return str(s).split(splitter)[i]
  else:
    return str(s).split(splitter)[i:j]


def process_summary(l, scans):
  i = None
  if l.startswith('qmass'):
    lhs, rhs = piece(l, "=", 0, 2)
    i = int(lhs[5:])
    mass = rhs
    scans.setdefault(i, {})['theory_mass'] = float(mass)
  if l.startswith('qexp'):
    lhs, rhs = piece(l, "=", 0, 2)
    i = int(lhs[4:])
    r, c = rhs.split(",")
    scans.setdefault(i, {})['exp_mass_charge'] = float(r)
    scans.setdefault(i, {})['charge0'] = c
  if l.startswith('qplughole'):
    lhs, rhs = l.strip().split("=")
    i = int(lhs[9:])
    scans.setdefault(i, {})['homology'] = float(rhs)
  if l.startswith('qmatch'):
    lhs, rhs = l.strip().split("=")
    i = int(lhs[6:])
    val = float(rhs)
    if val == 0:
      identity = 0.0
    else:
      identity = 10.0*math.log10(float(val))
    scans.setdefault(i, {})['identity'] = identity
  if i is not None: 
    if 'matches' not in scans[i]:
      scans[i]['matches'] = []


def process_matches(l, scans, max_peptide_rank):
  query, data_str = piece(l, "=", 0, 2)
  if data_str.strip() == "-1":
    return
  if 'subs' in query:
    return

  entry, match = piece(query, "_", 0, 2)
  if max_peptide_rank is not None:
    i_match = int(match[1:])
    if i_match > max_peptide_rank:
      return

  query_pieces = query.split("_")
  scan_id = int(entry[1:])
  if 'terms' in query:
    if scan_id not in scans:
      return
    terms_list = data_str.split(":")
    sites = scans[scan_id]['matches'][-1]['proteins']
    for terms, site in zip(terms_list, sites):
      n_term, c_term = terms.split(",")
      site['n_term'] = n_term
      site['c_term'] = c_term
  elif len(query_pieces) == 2:
    d = data_str.split(",")
    try:
      data = {
        'skip': int(d[0]),
        'mass': float(d[1]),
        'mass_diff': float(d[2]),
        'sequence': d[4],
        'score': float(d[7]),
        'proteins': []
      }
    except:
      raise ValueError("Error reading: " + l)
    protein_sites_str = piece(l, ";", 1)
    protein_sites = protein_sites_str.split(",")

    i_match = len(scans[scan_id]['matches'])
    score = data['score']
    for protein_site in protein_sites:
      p = protein_site.split(":")
      seqid = p[0][1:-1]
      protein = {
        'seqid': seqid,
        'i': int(p[2]),
        'j': int(p[3]),
      }
      data['proteins'].append(protein)
    scans[scan_id]['matches'].append(data)


def process_proteins(l, proteins):
  i = l.index(",")
  lhs = l[:i]
  rhs = l[i:]
  words = lhs.split("=")
  seqid = words[0][1:-1]
  mass = float(words[1])
  name = rhs[1:-1]
  proteins[seqid] = {
    'mass': mass,
    'description': name
  }


def process_query(l, scan_id, scans):
  lhs, rhs = l.split("=")
  if lhs in scans[scan_id]:
    raise ValueError(
        "%s already in scans[%d]" % (lhs, scan_id))
  if lhs == 'title':
    scans[scan_id][lhs] = urllib.unquote(rhs)
  else:
    scans[scan_id][lhs] = parse_string(rhs)


def read_mascot_dat(fname, max_peptide_rank=1):
  scans = {}
  proteins = {}
  boundary_id = None
  section = None
  process_line = None
  for l in open(fname, 'rU'):
    if not l.strip():
      continue
    if not boundary_id:
      if l.startswith("Content-Type") and 'boundary' in l:
        boundary_id = "--" + piece(l, "=", -1)
        continue
    if not boundary_id:
      continue
    if boundary_id and l.startswith(boundary_id):
      section = None
      process_line = None
      continue
    if l.startswith("Content-Type"):
      word = piece(l, None, -1)
      name = piece(word, "=", -1)
      section = name[1:-1]
      if section == "summary":
        process_line = \
            lambda l: process_summary(l, scans)
      if section == "peptides":
        process_line = \
            lambda l: process_matches(l, scans, max_peptide_rank)
      if section == "proteins":
        process_line = \
            lambda l: process_proteins(l, proteins)
      if "query" in section:
        scan_id = int(section[5:])
        process_line = \
            lambda l: process_query(l, scan_id, scans)
      continue
    if process_line:
      process_line(l[:-1])
  return scans, proteins


def split_mascot_ion_str(s):
    "Parses an Ion entry string into a dictionary"
    pairs = [piece.split(':') for piece in s.split(',')]
    return [[float(x),float(y)] for x,y in pairs]


def load_mascot_dat_to_proteins(proteins, i_source, mascot_dat):
  peptide_by_match_id = {}
  for seqid in proteins:
    protein = proteins[seqid]
    for source in protein['sources']:
      for peptide in source['peptides']:
        if 'identity' in peptide['attr']:
          match_id = "%.2f%.2f%.2f%s" % \
              (peptide['attr']['score'],
               peptide['attr']['identity'],
               peptide['attr']['homology'],
               peptide['sequence'])
          peptide_by_match_id[match_id] = peptide
  scans, mascot_proteins = read_mascot_dat(mascot_dat)
  n_match = 0
  for scan in scans.values():
    for match in scan['matches']:
      match_id = "%.2f%.2f%.2f%s" % \
          (match['score'],
           scan['identity'],
           scan['homology'],
           match['sequence'])
      if match_id in peptide_by_match_id:
        n_match += 1
        peptide = peptide_by_match_id[match_id]
        peptide['spectrum'] = split_mascot_ion_str(scan['Ions1'])
  logger.info('%s: matched %d pepXML to %d mascot PSM' % \
      (mascot_dat, len(peptide_by_match_id), n_match))


if __name__ == '__main__':
  scans, proteins = read_mascot_dat('../example/mascot/F022045.dat')
  save_data_dict(scans, '../example/mascot/scans.dump')
  save_data_dict(proteins,'../example/mascot/proteins.dump')
