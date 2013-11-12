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
import logging

import uniprot
import tpp
import parse
import proteins as peptagram_proteins


logger = logging.getLogger('protxml2csv')


usage = """
Reads a TPP Protein XML file and outputs it into a convenient CSV format.
For better protein name description, interrogates the Uniprot website to
extract more useful meta-data.

Usage: %s protein_xml [output_csv]
""" % os.path.basename(__file__)


def clean_seqid(seqid):
  if '|' in seqid:
    return seqid.split('|')[1]
  else:
    return seqid


def get_all_seqids(proteins):
  results = []
  for seqid in proteins:
    results.append(seqid)
    results.extend(proteins[seqid]['attr']['other_seqids'])
  return results


def is_url_connected(url):
  "Tests if a connection can be made to url"
  try:
    response=urllib2.urlopen(url, timeout=1)
    return True
  except urllib2.URLError as err: 
    pass
  return False


def add_uniprot_data(proteins, cache_basename=None):
  """
  Processes the data from an PROTXML file, reads the
  seqids, and attempts to mapt to a UNIPROT ID
  and then, fetch the metadata of that protein from
  the uniprot.org website: organism, gene, description
  etc.
  """
  seqids = get_all_seqids(proteins)
  if is_url_connected('http://uniprot.org'):
    uniprot_dict = uniprot.get_metadata_with_some_seqid_conversions(seqids, cache_basename)
  else:
    print "Can't connect to www.uniprot.org, won't use uniprot metatdata"
    uniprot_dict = {}
  for seqid in proteins:
    protein = proteins[seqid]
    if 'other_seqids' not in protein['attr']:
      new_seqid = seqid
    else:
      names = [seqid] + protein['attr']['other_seqids']
      new_seqid = uniprot.sort_seqids_by_uniprot(names, uniprot_dict)[0]
      if new_seqid != seqid:
        print "Protein is better represented with %s than %s" % \
            (uniprot.get_naked_seqid(new_seqid), 
             uniprot.get_naked_seqid(seqid))
        protein['attr']['seqid'] = new_seqid
        protein['attr']['other_seqids'] = names[1:]
        proteins[new_seqid] = protein
        del proteins[seqid]
    if new_seqid not in uniprot_dict:
      print "No uniprot metadata for seqid %s" % \
          (uniprot.get_naked_seqid(new_seqid))
      continue
    uniprot_entry = uniprot_dict[new_seqid]
    protein['attr']['uniprot_id'] = uniprot_entry['id']
    protein['attr']['seqid'] = new_seqid
    protein['attr']['link'] = \
        '=HYPERLINK("http://uniprot.org/uniprot/%s")' % \
            uniprot_entry['id']
    if 'gene' in uniprot_entry:
      protein['attr']['gene'] = uniprot_entry['gene']
    if 'organism' in uniprot_entry:
      protein['attr']['organism'] = uniprot_entry['organism']
    protein['description'] = '; '.join(uniprot_entry['descriptions'])
    if 'Uncharacterized protein' in protein['description']: 
      if 'comment' in uniprot_entry:
        if '-!- SIMILARITY' in uniprot_entry['comment']:
          found_similarity = False
          protein['description'] = ''
          for line in uniprot_entry['comment'].splitlines():
            if found_similarity and '-!-' in line:
              break
            if '-!- SIMILARITY' in line:
              found_similarity = True
            protein['description'] += line.strip()
    protein['sequence'] = uniprot_entry['sequence']
    protein['attr']['length'] = len(protein['sequence'])


def write_proteins_to_csv(proteins, proteins_csv):
  headings = [
    'group_id',
    'sibling',
    'seqid',
    'uniprot_id',
    'link',
    'description',
    'gene',
    'length',
    'percent_coverage',
    'probability',
    'n_peptide',
    'n_spectrum',
    'percent_spectra',
    'n_unique_peptide',
    'n_unique_spectrum',
    'percent_unique_spectra',
    'organism',
    'other_seqids'
    ]

  rows = []
  for seqid in proteins:
    protein = proteins[seqid]
    row = []
    for key in headings:
      if key in protein['attr']:
        if key == 'other_seqids':
          val = ','.join(protein['attr'][key])
        else:
          val = protein['attr'][key]
      elif key in protein:
        val = protein[key]
      else:
        val = ''
      row.append(val)
    rows.append(row)
  rows.sort(key=lambda r:r[0])
  rows.insert(0, headings)

  print "-"*60
  print "WRITING:", os.path.abspath(proteins_csv)
  csv_writer = csv.writer(open(proteins_csv, 'wb'))
  for row in rows:
    csv_writer.writerow(row)


def convert_to_csv(params):
  proteins, sources = tpp.get_proteins_and_sources(
      params['protxml'], 
      params['pepxmls'], 
      [params['peptide_error_cutoff']],
      params['protein_error_cutoff'])

  clean_seqid = None
  if 'clean_seqid_fn' in params:
    clean_seqid = params['clean_seqid_fn']
  peptagram_proteins.change_seqids_in_proteins(proteins, clean_seqid)

  cache_basename = params['protxml'].replace('.prot.xml', '.uniprot')
  add_uniprot_data(proteins, cache_basename)

  proteins_csv = params['csv']
  if not proteins_csv:
    proteins_csv = params['protxml'] + '.csv'
  write_proteins_to_csv(proteins, proteins_csv)


if __name__ == "__main__":
  params = {
    'protxml': '../data/kl4/TPP/LNCAP_Exp1/interactLNCap-KLK4-1-10.prot.xml',
    'pepxmls': ['../data/kl4/TPP/LNCAP_Exp1/interactlncap-klk4-1-10.pep.xml'],
    'csv': '../data/kl4/TPP/LNCAP_Exp1/interactLNCap-KLK4-1-10.csv',
    'protein_error_cutoff': 0.01, # or None
    'peptide_probability_cutoff': 0.5, # or None
    'peptide_error_cutoff': None, # or None
    'clean_seqid_fn': clean_seqid,
  }
  convert_to_csv(params)




