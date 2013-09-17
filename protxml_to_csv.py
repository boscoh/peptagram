# -*- coding: utf-8 -*-
# (c) 2013 Bosco Ho
# from __future__ import print_function
import os
import sys
import csv
import json
from pprint import pprint
import re
import urllib2

import uniprot
import peptagram.tpp
import peptagram.parse


usage = """
Reads a TPP Protein XML file and outputs it into a convenient CSV format.
For better protein name description, interrogates the Uniprot website to
extract more useful meta-data.

Usage: %s protein_xml [output_csv]
""" % os.path.basename(__file__)


def get_all_seqids(protxml_groups):
  """
  Returns all seqids, including alternative seqids from
  the data structure parsed from an PROTXML file.
  """
  results = []
  for group_id, protxml_group in protxml_groups.items():
    for protxml_protein in protxml_group['proteins']:
      results.append(protxml_protein['protein_name'])
      results.extend(protxml_protein['other_seqids'])
  return results


def is_url_connected(url):
  "Tests if a connection can be made to url"
  try:
    response=urllib2.urlopen(url, timeout=1)
    return True
  except urllib2.URLError as err: 
    pass
  return False


def add_uniprot_data(protxml_groups, cache_file=None):
  """
  Processes the data from an PROTXML file, reads the
  seqids, and attempts to mapt to a UNIPROT ID
  and then, fetch the metadata of that protein from
  the uniprot.org website: organism, gene, description
  etc.
  """
  seqids = get_all_seqids(protxml_groups)
  if is_url_connected('http://uniprot.org'):
    uniprot_dict = uniprot.get_metadata_with_some_seqid_conversions(seqids, cache_file)
  else:
    print "Can't connect to www.uniprot.org, won't use uniprot metatdata"
    uniprot_dict = {}
  for group_id, protxml_group in protxml_groups.items():
    for protein in protxml_group['proteins']:
      protein['id'] = ''
      protein['acc'] = protein['protein_name']
      names = [protein['protein_name']] + protein['other_seqids']
      new_seqid = uniprot.sort_seqids_by_uniprot(names, uniprot_dict)[0]
      if new_seqid != protein['protein_name']:
        print "Protein group %s%s is better represented with %s than %s" % \
            (group_id, 
             protein['group_sibling_id'],
             uniprot.get_naked_seqid(new_seqid), 
             uniprot.get_naked_seqid(protein['protein_name']))
        protein['protein_name'] = new_seqid
        protein['other_seqids'] = names[1:]
        protein['acc'] = new_seqid
      protein['other_seqids'] = ';'.join(protein['other_seqids'])
      if new_seqid not in uniprot_dict:
        print "No uniprot metadata for protein group %s%s seqid %s" % \
            (group_id, 
             protein['group_sibling_id'],
             uniprot.get_naked_seqid(new_seqid))
        continue
      protein['link'] = ''
      uniprot_entry = uniprot_dict[new_seqid]
      protein['id'] = uniprot_entry['id']
      protein['acc'] = uniprot_entry['accs'][0]
      protein['link'] = \
          '=HYPERLINK("http://uniprot.org/uniprot/%s")' % \
              uniprot_dict[new_seqid]['id']
      if 'gene' in uniprot_entry:
        protein['gene'] = uniprot_entry['gene']
      if 'organism' in uniprot_entry:
        protein['organism'] = uniprot_entry['organism']
      protein['description'] = '; '.join(uniprot_entry['descriptions'])
      if 'length' in uniprot_entry:
        protein['length'] = uniprot_entry['length']


def convert_to_csv(params):
  protxml_fname = params['protxml']
  proteins_csv = params['csv']
  if not proteins_csv:
    proteins_csv = protxml_fname + '.csv'
  is_skip_no_unique = params['is_skip_no_unique']  
  protein_probability_cutoff = params['protein_probability_cutoff']
  protein_error_cutoff = params['protein_error_cutoff']

  print "LOADING:", protxml_fname
  protxml_groups, protein_probs = peptagram.tpp.read_protxml(protxml_fname)

  protein_by_seqid = {}
  for group_id, protxml_group in protxml_groups.items():
    for protein in protxml_group['proteins']:
       seqids = [protein['protein_name']] + protein['other_seqids']
       for seqid in seqids:
         protein_by_seqid[seqid] = protein

  if 'pepxmls' in params:
    n_scan = 0
    for pepxml in params['pepxmls']:
      print pepxml
      pepxml_reader = peptagram.tpp.PepxmlReader(
          pepxml, error_cutoff=params['peptide_error_cutoff'])
      for scan in pepxml_reader.iter():
        for match in scan['matches']:
          seqid = match['protein']
          if seqid not in protein_by_seqid:
            pprint('Warning {}'.format(seqid))
            continue
          protein = protein_by_seqid[seqid]
          if 'n_scan' not in protein:
            protein['n_scan'] = 0
            protein['scanned_peptides'] = set()
          protein['scanned_peptides'].add(match['peptide'])
          protein['n_scan'] += 1
        n_scan += 1
    for group_id, protxml_group in protxml_groups.items():
      for protein in protxml_group['proteins']:
        if 'n_scan' not in protein:
          print protein['protein_name'], 'has no scans'
          continue
        percent = 100.0*protein['n_scan']/n_scan
        percent = float('%.2f' % percent)
        protein['percent_spectra'] = percent
        n_peptide = 0
        n_peptide_unique = 0
        for peptide in protein['peptides']:
          if peptide['peptide_sequence'] in protein['scanned_peptides']:
            n_peptide += 1
            if peptide['is_nondegenerate_evidence'] == 'Y':
              n_peptide_unique += 1
        protein['n_peptide'] = n_peptide
        protein['n_unique_peptide'] = n_peptide_unique

  prob_cutoff = peptagram.tpp.error_to_probability(protein_probs, protein_error_cutoff)
  print('protein probability cutoff', prob_cutoff)
  for group_id in protxml_groups.keys():
    protxml_group = protxml_groups[group_id]
    proteins = protxml_group['proteins']
    for i in reversed(range(len(proteins))):
      protein = proteins[i]
      if protein['probability'] < prob_cutoff:
        del proteins[i]
      elif 'n_scan' not in protein or protein['n_scan'] == 0:
        del proteins[i]
      elif 'n_peptide' not in protein or protein['n_peptide'] == 0:
        del proteins[i]
    if len(proteins) == 0:
      del protxml_groups[group_id]

  cache_file = proteins_csv.replace('.csv', '.uniprot')
  add_uniprot_data(protxml_groups, cache_file)

  title_key_pairs = [
      ('group', 'group_number'),
      ('sibling', 'group_sibling_id'),
      ('seqid', 'acc'),
      ('uniprot_id', 'id'),
      ('url', 'link'),
      ('description', 'description'),
      ('gene', 'gene'),
      ('length', 'length'),
      ('percent_coverage', 'percent_coverage'),
      ('probability', 'probability'),
      ('num_peptide', 'n_peptide'),
      ('num_nondegenerate_peptide', 'n_unique_peptide'),
      ('num_spectra', 'n_scan'),
      ('percent_spectra', 'percent_spectra'),
      ('organism', 'organism'),
      ('other_seqids', 'other_seqids')
      ]

  headings = [title for title, key in title_key_pairs]
  rows = []
  for group_id, group in protxml_groups.items():
    for protein in group['proteins']:
      row = []
      for title, key in title_key_pairs:
        row.append(protein[key] if key in protein else "")
      rows.append(row)
  rows.sort(key=lambda r:r[0])
  rows.insert(0, headings)

  print "-"*60
  print "WRITING:", os.path.abspath(proteins_csv)
  csv_writer = csv.writer(open(proteins_csv, 'wb'))
  for row in rows:
    csv_writer.writerow(row)


if __name__ == "__main__":
  params = {
    'protxml': 'data/kl4/TPP/protxml/interactWMPY-1.prot.xml',
    'pepxmls': ['data/kl4/TPP/protxml/interactwmpy-1.pep.xml'],
    'csv': 'data/kl4/TPP/protxml/interactwmpy-1.csv',
    'protein_probability_cutoff': None, # or None
    'protein_error_cutoff': 0.01, # or None
    'peptide_probability_cutoff': 0.5, # or None
    'peptide_error_cutoff': None, # or None
    'is_skip_no_unique': False,
  }
  convert_to_csv(params)




