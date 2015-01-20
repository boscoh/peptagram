# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import math
import urllib
import xml.etree.ElementTree as etree

import logging
logger = logging.getLogger('mascot')

from parse import parse_string, save_data_dict
import proteins as parse_proteins


"""
Parser for Mascot dat files 

Some useful formulae:
 - identity_score = -10*Math.log10(1/qmatch)
 - qplughole = homology_score
 - significant if identity_score > homology_score

Main API entry:

  get_proteins(
      mascot_dat, 
      great_score=80, 
      cutoff_score=0)

  returns a dictionary that organizes peptide-spectrum-matches
  around proteins, and a list of sources.

"""


def splitter(s, delimiter, i, j=None):
  if j is None:
    return str(s).split(delimiter)[i]
  else:
    return str(s).split(delimiter)[i:j]


class MascotReader():
  def __init__(self, mascot_dat, max_peptide_rank=1):
    self.mascot_dat = mascot_dat
    self.scans = {}
    self.proteins = {}
    self.masses = {}
    self.max_peptide_rank = max_peptide_rank
    self.unimod_lines = []
    self.unimod = None
    self.scan_id = None

    self.section = None
    self.boundary_id = None

    self.process_line = None
    self.read_mascot_dat()

  def read_mascot_dat(self):
    for l in open(self.mascot_dat):
      if not l.strip():
        continue

      if not self.boundary_id:
        if l.startswith("Content-Type") and 'boundary' in l:
          self.boundary_id = "--" + splitter(l, "=", -1)
          continue

      if not self.boundary_id:
        continue

      if self.boundary_id and l.startswith(self.boundary_id):
        self.section = None
        self.process_line = None
        continue

      if l.startswith("Content-Type"):
        word = splitter(l, None, -1)
        name = splitter(word, "=", -1)
        self.section = name[1:-1]
        if self.section == "masses":
          self.process_line = self.process_masses
        if self.section == "unimod":
          self.process_line = self.process_unimod
        if self.section == "summary":
          self.process_line = self.process_summary
        if self.section in ['matches', 'peptides']:
          self.process_line = self.process_matches
        if "query" in self.section:
          self.scan_id = int(self.section[5:])
          self.process_line = self.process_query
        if self.section == "proteins":
          self.process_line = self.process_proteins
        continue

      if self.process_line:
        self.process_line(l[:-1])

  def process_summary(self, l):
    i = None
    if l.startswith('qmass'):
      lhs, rhs = splitter(l, "=", 0, 2)
      i = int(lhs[5:])
      mass = rhs
      self.scans.setdefault(i, {})['theory_mass'] = float(mass)
    if l.startswith('qexp'):
      lhs, rhs = splitter(l, "=", 0, 2)
      i = int(lhs[4:])
      r, c = rhs.split(",")
      self.scans.setdefault(i, {})['exp_mass_charge'] = float(r)
      self.scans.setdefault(i, {})['charge0'] = c
    if l.startswith('qplughole'):
      lhs, rhs = l.strip().split("=")
      i = int(lhs[9:])
      self.scans.setdefault(i, {})['homology'] = float(rhs)
    if l.startswith('qmatch'):
      lhs, rhs = l.strip().split("=")
      i = int(lhs[6:])
      val = float(rhs)
      if val == 0:
        identity = 0.0
      else:
        identity = 10.0*math.log10(float(val))
      self.scans.setdefault(i, {})['identity'] = identity
    if i is not None: 
      if 'matches' not in self.scans[i]:
        self.scans[i]['matches'] = []

  def process_matches(self, l):
    query, data_str = splitter(l, "=", 0, 2)
    if data_str.strip() == "-1":
      return
    if 'subs' in query:
      return

    entry, match = splitter(query, "_", 0, 2)
    if self.max_peptide_rank is not None:
      i_match = int(match[1:])
      if i_match > self.max_peptide_rank:
        return

    query_pieces = query.split("_")
    scan_id = int(entry[1:])
    if 'terms' in query:
      if scan_id not in self.scans:
        return
      terms_list = data_str.split(":")
      sites = self.scans[scan_id]['matches'][-1]['proteins']
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
          'mod_mask_str': d[6],
          'score': float(d[7]),
          'proteins': []
        }
      except:
        raise ValueError("Error reading: " + l)
      protein_sites_str = splitter(l, ";", 1)
      protein_sites = protein_sites_str.split(",")

      i_match = len(self.scans[scan_id]['matches'])
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
      self.scans[scan_id]['matches'].append(data)

  def process_proteins(self, l):
    i = l.index(",")
    lhs = l[:i]
    rhs = l[i:]
    words = lhs.split("=")
    seqid = words[0][1:-1]
    mass = float(words[1])
    name = rhs[1:-1]
    self.proteins[seqid] = {
      'mass': mass,
      'description': name
    }

  def process_unimod(self, l):
    if l.startswith("<?xml"):
      self.unimod_lines = []
    self.unimod_lines.append(l)
    if not l.startswith("</umod:unimod>"):
      return
    unimod_text = '\n'.join(self.unimod_lines)
    self.unimod = etree.fromstring(unimod_text)
    
  def process_masses(self, l):
    if not l:
      return
    l, r = l.split('=')
    if l.startswith('delta'):
      num = l.replace('delta', '')
      delta_mass, name = r.split(',')
      self.masses[num] = {
        'delta_mass': float(delta_mass),
        'type': name
      }
    elif len(l) == 1:
      self.masses[l] = float(r)

  def process_query(self, l):
    lhs, rhs = l.split("=")
    if lhs in self.scans[self.scan_id]:
      raise ValueError(
          "%s already in scans[%d]" % (lhs, self.scan_id))
    if lhs == 'title':
      self.scans[self.scan_id]['title'] = urllib.unquote(rhs)
    elif lhs == 'Ions1':
      pairs = [piece.split(':') for piece in rhs.split(',')]
      ions = [[float(x),float(y)] for x,y in pairs]
      self.scans[self.scan_id]['Ions1'] = ions
    else:
      self.scans[self.scan_id][lhs] = parse_string(rhs)


def get_proteins(mascot_dat, great_score=80, cutoff_score=0):
  mascot_reader = MascotReader(mascot_dat)
  scans = mascot_reader.scans
  mascot_proteins = mascot_reader.proteins
  masses = mascot_reader.masses

  proteins = {}
  for scan_id, scan in scans.items():
    for mascot_match in scan['matches']:
      peptide_sequence = mascot_match['sequence']
      match = parse_proteins.new_match(peptide_sequence)
      match['attr']['missed_cleavages'] = mascot_match['skip']
      match['ionscore'] = mascot_match['score']
      match['spectrum'] = scan['Ions1']
      score = match['ionscore']
      intensity = \
          parse_proteins.calc_intensity(
              score, great_score, cutoff_score)
      match['intensity'] = intensity
      for key in scan:
        if key not in ["Ions1", "matches"]:
          match['attr'][key] = scan[key]
      n = len(peptide_sequence)
      match['modifications'] = []
      for i_peptide in range(-1, n+1):
        num = mascot_match['mod_mask_str'][i_peptide+1]
        if num in masses:
          delta_mass = masses[num]['delta_mass']
          mass = masses[peptide_sequence[i_peptide]] + delta_mass
          modification = {
            'i': i_peptide,
            'mass': mass,
            'type': masses[num]['type']
          }
          match['modifications'].append(modification)
      for match_protein in mascot_match['proteins']:
        seqid = match_protein['seqid']
        if seqid not in proteins:
          proteins[seqid] = parse_proteins.new_protein(seqid)
        protein = proteins[seqid]
        protein['sources'][0]['matches'].append(match)

  return proteins


# def load_mascot_dat_to_proteins(proteins, i_source, mascot_dat):
#   peptide_by_match_id = {}
#   for seqid in proteins:
#     protein = proteins[seqid]
#     for source in protein['sources']:
#       for peptide in source['matches']:
#         if 'identity' in peptide['attr']:
#           match_id = "%.2f%.2f%.2f%s" % \
#               (peptide['attr']['score'],
#                peptide['attr']['identity'],
#                peptide['attr']['homology'],
#                peptide['sequence'])
#           peptide_by_match_id[match_id] = peptide
#   scans, mascot_proteins, masses = read_mascot_dat(mascot_dat)
#   n_match = 0
#   for scan in scans.values():
#     for match in scan['matches']:
#       match_id = "%.2f%.2f%.2f%s" % \
#           (match['score'],
#            scan['identity'],
#            scan['homology'],
#            match['sequence'])
#       if match_id in peptide_by_match_id:
#         n_match += 1
#         peptide = peptide_by_match_id[match_id]
#         peptide['spectrum'] = split_mascot_ion_str(scan['Ions1'])
#   logger.info('%s: matched %d pepXML to %d mascot PSM' % \
#       (mascot_dat, len(peptide_by_match_id), n_match))
