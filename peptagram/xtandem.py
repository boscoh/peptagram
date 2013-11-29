# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint
import xml.etree.ElementTree as etree
import json
import os
import re
import logging

import parse
import proteins as proteins_module


"""
Parser for X!Tandem XML output.
"""


logger = logging.getLogger('xtandem')


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
      'description': description
    }
    domain_elem = peptide_elem.find('domain')
    match.update(parse.parse_attrib(domain_elem))
    scan['matches'].append(match)

  elem = top_elem.find('group[@label="fragment ion mass spectrum"]')

  note = elem.find('note')
  if note:
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


def load_xtandem_into_proteins(proteins, xtandem_fname, i_source, n_peak=50):
  proteins_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    proteins_by_seqid[seqid] = protein
    for other_seqid in protein['attr']['other_seqids']:
      proteins_by_seqid[other_seqid] = protein

  fastas = {}

  for scan in read_xtandem(xtandem_fname):
    scan_id = scan['id']
    is_match = False
    for match in scan['matches']:
      seqid = match['seqid']
      if seqid not in proteins_by_seqid:
        continue
      protein = proteins_by_seqid[seqid]
      if seqid not in fastas:
        fastas[seqid] = {
          'sequence': match['sequence'],
          'description': match['description'],
        }
      source = protein['sources'][i_source]    
      for peptide in source['peptides']:
        if scan_id != peptide['attr']['scan_id']:
          continue
        is_match = True
        x_vals = map(float, scan['masses'].split())
        y_vals = map(float, scan['intensities'].split())
        ions = [(x, y) for x, y in zip(x_vals, y_vals)]
        ions.sort(key=lambda i:-i[1])
        peptide['spectrum'] = ions[:n_peak]

  proteins_module.load_fastas_into_proteins(proteins, fastas)
  for seqid in proteins.keys():
    protein = proteins[seqid]
    if 'sequence' not in protein:
      logger.debug("Protein {} not found in x!tandem".format(seqid))
      del proteins[seqid]
      continue
    n_peptide = sum([len(source['peptides']) for source in protein['sources']])
    if n_peptide == 0:
      del proteins[seqid]
      logger.debug("No peptide-spectra matches found in {}".format(seqid))
      continue

  # proteins_module.calculate_peptide_positions(proteins)



if __name__ == "__main__":
  scans, fastas = read('../example/xtandem/Seq23282_E1O1.tandem')
  parse.save_data_dict(scans, '../example/xtandem/scans.dump')
  parse.save_data_dict(fastas, '../example/xtandem/fastas.dump')





