# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint
import xml.etree.ElementTree as etree
import json
import os
import re
import logging

import parse
import mass



"""
Parser for X!Tandem XML output.
"""


logger = logging.getLogger('xtandem')



def strip_whitespace(txt):
  result = ""
  for line in txt.splitlines():
    result += ''.join(line.split())
  return result


def parse_scan(group, nsmap):
  scan = {}
  fastas = {}

  scan.update(parse.parse_attrib(group))

  scan['matches'] = []
  for protein in group.findall('protein'):
    match = {}

    words = protein.attrib['label'].split()
    seqid = words[0]
    description = ' '.join(words[1:])
    peptide_elem = protein.find('peptide')

    match['seqid'] = seqid

    if seqid not in fastas:
      sequence = strip_whitespace(peptide_elem.text)
      fastas[seqid] = {
        'description': description,
        'sequence': sequence,
      }

    domain_elem = peptide_elem.find('domain')
    match.update(parse.parse_attrib(domain_elem))
    scan['matches'].append(match)

  elem = group.find('group[@label="fragment ion mass spectrum"]')

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

  return scan, fastas


def read(xtandem_xml):
  fastas = {}
  scans = []
  nsmap = {}
  for event, elem in etree.iterparse(xtandem_xml, events=('end', 'start-ns')):
    if event == 'start-ns':
      nsmap.update({elem})
    if event == 'end':
      if elem.tag == 'group' and elem.attrib['type'] == 'model':
        scan, scan_fastas = parse_scan(elem, nsmap)
        scans.append(scan)
        fastas.update(scan_fastas)
        elem.clear()
  return scans, fastas


def merge_peaks(unmatched_peaks, matched_peaks):
  peaks = []
  for m, i in unmatched_peaks:
    peaks.append([m, i, ''])
  peaks.extend(matched_peaks)
  return peaks


def load_scans_into_proteins(proteins, xtandem_scans, i_source, n_peak=50):
  """
  Only takes top n_peak from xtandem.
  """
  scans = { int(scan['id']):scan for scan in xtandem_scans }
  for seqid in proteins.keys():
    source = proteins[seqid]['sources'][i_source]
    scan_ids = []
    for peptide in source['peptides']:
      scan_id = peptide['attr']['scan_id']
      sequence = peptide['sequence']
      modifications = peptide['attr']['modifications']
      if scan_id not in scans:
        logger.warning("Couldn't find xtadnem entry for scan {} in pepxml".format(scan_id))
        continue
      scan = scans[scan_id]
      x_vals = map(float, scan['masses'].split())
      y_vals = map(float, scan['intensities'].split())
      ions = [(x, y) for x, y in zip(x_vals, y_vals)]
      ions.sort(key=lambda i:-i[1])
      peptide['spectrum'] = ions[:n_peak]



if __name__ == "__main__":
  scans, fastas = read('../example/xtandem/Seq23282_E1O1.tandem')
  parse.save_data_dict(scans, '../example/xtandem/scans.dump')
  parse.save_data_dict(fastas, '../example/xtandem/fastas.dump')





