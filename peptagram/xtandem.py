# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import re
import math
import xml.etree.ElementTree as etree

import logging
logger = logging.getLogger('xtandem')

import parse
import proteins as proteins_module
import peptidemass


"""
Parser for X!Tandem XML mass-spec search results.

Main API entry:

  get_proteins(
    xtandem_fname, 
    n_peak=50, 
    good_expect=1E-8,
    cutoff_expect=1E-2, 
    excluded_seqids=[], 
    include_seqids=[])

  returns a dictionary that organizes peptide-spectrum-matches
  around proteins.

"""


def strip_whitespace(txt):
  result = ""
  for line in txt.splitlines():
    result += ''.join(line.split())
  return result


def parse_scan(top_elem, nsmap):
  scan = {}
  scan.update(parse.parse_attrib(top_elem))

  scan['matches'] = []
  for protein in top_elem.findall('protein'):
    words = protein.attrib['label'].split()
    seqid = words[0]
    description = ' '.join(words[1:])
    peptide_elem = protein.find('peptide')
    sequence = strip_whitespace(peptide_elem.text)
    match = {  
      'seqid': seqid,
      'sequence': sequence,
      'description': description,
      'modifications': []
    }
    domain_elem = peptide_elem.find('domain')
    for mod_elem in domain_elem.findall('aa'):
      match['modifications'].append(parse.parse_attrib(mod_elem))
    match.update(parse.parse_attrib(domain_elem))
    scan['matches'].append(match)

  elem = top_elem.find('group[@label="fragment ion mass spectrum"]')

  note = elem.find('note')
  if note is not None:
    scan['Description'] = elem.find('note').text.strip()
  else:
    scan['Description'] = ''

  masses_elem = elem.find(
      'GAML:trace/GAML:Xdata/GAML:values', namespaces=nsmap)
  scan['masses'] = masses_elem.text.strip()

  intensities_elem = elem.find(
      'GAML:trace/GAML:Ydata/GAML:values', namespaces=nsmap)
  scan['intensities'] = intensities_elem.text.strip()

  charge_elem = elem.find(
      'GAML:trace/GAML:attribute[@type="charge"]', namespaces=nsmap)
  scan['charge'] = charge_elem.text.strip()

  mass_elem = elem.find(
      'GAML:trace/GAML:attribute[@type="M+H"]', namespaces=nsmap)
  scan['mass'] = mass_elem.text.strip()

  return scan


def read_xtandem(xtandem_xml):
  fastas = {}
  scans = []
  nsmap = {}
  for event, elem in etree.iterparse(xtandem_xml, events=('end', 'start-ns')):
    if event == 'start-ns':
      nsmap.update({elem})
    if event == 'end':
      if elem.tag == 'group' and elem.attrib['type'] == 'model':
        scan = parse_scan(elem, nsmap)
        yield scan
        elem.clear()


def get_proteins(
    xtandem_fname, 
    n_peak=50, 
    good_expect=1E-8,
    cutoff_expect=1E-2):
  proteins = {}
  i_source = 0
  print_scan = True
  for scan in read_xtandem(xtandem_fname):
    scan_id = scan['id']
    x_vals = map(float, scan['masses'].split())
    y_vals = map(float, scan['intensities'].split())
    ions = [(x, y) for x, y in zip(x_vals, y_vals)]
    ions.sort(key=lambda i:-i[1])

    for xtandem_match in scan['matches']:

      expect = xtandem_match['expect']
      if cutoff_expect < expect:
        continue

      intensity = proteins_module.calc_minus_log_intensity(
        expect, good_expect, cutoff_expect)

      seqid = xtandem_match['seqid']

      if seqid not in proteins:
        protein = proteins_module.new_protein(seqid)
        protein.update({
          'sequence': xtandem_match['sequence'],
          'description': xtandem_match['description'],
        })
        proteins[seqid] = protein

      protein = proteins[seqid]
      source = protein['sources'][i_source]    

      match = {
        'sequence': xtandem_match['seq'],
        'intensity': intensity,
        'modifications': [],
        'mask': mask,
        'spectrum': ions[:n_peak],
        'attr': {
          'scan_id': scan['id'],
          'charge': scan['charge'],
          'expect': expect,
          'missed_cleavages': xtandem_match['missed_cleavages'],
          'mass': scan['mass'],
          'source': parse.basename(xtandem_fname),
        }
      }
      if xtandem_match['modifications']:
        for mod in xtandem_match['modifications']:
          i_mod_in_full_seq = int(mod['at'])-1
          full_seq = xtandem_match['sequence']
          i_pep_seq = int(xtandem_match['start'])-1
          aa = mod['type']
          if aa in peptidemass.aa_monoisotopic_mass:
            mass = peptidemass.aa_monoisotopic_mass[aa]
          else:
            mass = 0.0
          match['modifications'].append({
            'i': i_mod_in_full_seq - i_pep_seq,
            'mass': mod['modified'] + mass,
          })

      source['matches'].append(match)
      
  proteins_module.calculate_peptide_positions(proteins)

  return proteins




