 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import os
import json
import math

import numpy 

import parse
import proteins as parse_proteins

"""
Parser for Maxquant summary text files
"""


def read(in_dir):
  peptides = parse.read_tsv(os.path.join(in_dir, 'peptides.txt'))
  scans = parse.read_tsv(os.path.join(in_dir, 'msms.txt'))
  protein_groups = parse.read_tsv(os.path.join(in_dir, 'proteinGroups.txt'))
  evidence = parse.read_tsv(os.path.join(in_dir, 'evidence.txt'))
  return peptides, scans, protein_groups, evidence


def get_labeled_spectrum(scan):
  masses = parse.splitter(scan['masses'], float)
  intensities = parse.splitter(scan['intensities'], float)
  labels = parse.splitter(scan['matches'])
  return [x for x in zip(masses, intensities, labels)]


def transfer_attrs(source_dict, target_dict, parse_list):
  for key, convert_fn in parse_list:
    if key in source_dict:
      target_dict[key] = convert_fn(source_dict[key])


def float_or_none(float_str):
  if float_str == "":
    return None
  return float(float_str)


evidence_parse_list = [
  ('intensity', float_or_none),
  ('experiment', lambda s: s.lower()),
  ('ratio h/l', float_or_none),
  ('ratio h/l normalized', float_or_none),
]

peptide_parse_list = [
  ('ratio h/l variability [%]', float_or_none),
]

protein_parse_list = [
  ('ratio h/l', float),
  ('ratio h/l normalized', float),
  ('ratio h/l variability [%]', float),
]

scan_parse_list = [
  ('scan number', int),
  ('m/z', float),
  ('labeling state', int),
  ('retention time', float),
  ('pep', float),
]


def change_key(data, old_key, new_key):
  if not old_key in data or old_key == new_key:
    return
  data[new_key] = data[old_key]
  del data[old_key]


def get_proteins_and_sources(in_dir, is_leu_ile_isomeric=False,):
  peptide_list, scan_list, protein_group_list, evidence_list = \
        read(in_dir)

  peptides = { int(p['id']):p for p in peptide_list }
  scans = { int(s['id']):s for s in scan_list }
  protein_groups = { int(p['id']):p for p in protein_group_list }
  evidence_dict = { int(e['id']):e for e in evidence_list }

  parse.save_data_dict(peptides, in_dir + '/peptides.dump')
  parse.save_data_dict(scans, in_dir + '/scans.dump')
  parse.save_data_dict(protein_groups, in_dir + '/protein_groups.dump')
  parse.save_data_dict(evidence_dict, in_dir + '/evidence.dump')

  sources_set = set(e['raw file'] for e in evidence_dict.values())
  sources = [str(s) for s in sorted(sources_set)]
  i_sources = {source:k for k, source in enumerate(sources)}

  proteins = {}
  protein_by_group_id = {}
  for group_id, protein_group in protein_groups.items():
    protein = {
      'description': '',
      'attr': { 
        'group_id': group_id,
        'other_seqids': [],
      },
      'sources': [{ 'peptides': [] } for k in range(len(i_sources))],
    }
    transfer_attrs(protein_group, protein['attr'], protein_parse_list)

    seqids = parse.splitter(protein_group['protein ids'])
    proteins[seqids[0]] = protein
    protein['attr']['seqid'] = seqids[0]
    protein['attr']['other_seqids'] = seqids[1:]
    protein_by_group_id[group_id] = protein

  print("Matching sequences and scan in proteins")
  i_scan = 0
  n_scan = len(scans)
  for scan_id, scan in scans.items():
    i_scan += 1
    if i_scan % 5000 == 0:
      print("{}/{} scans processed".format(i_scan, n_scan))
    evidence_id = int(scan['evidence id'])
    evidence = evidence_dict[evidence_id]

    peptide_id = int(scan['peptide id'])
    peptide = peptides[peptide_id]
    for group_id in parse.splitter(str(scan['protein group ids'])):
      new_peptide = {
        'sequence': scan['sequence'],
        'spectrum': get_labeled_spectrum(scan),
        'attr' : {
          'modifications': [],
          'mq_scan_id': scan_id,
          'is_unique': peptide['unique (groups)'] == 'yes',
        }
      }
      transfer_attrs(scan, new_peptide['attr'], scan_parse_list)
      transfer_attrs(evidence, new_peptide['attr'], evidence_parse_list)
      transfer_attrs(peptide, new_peptide['attr'], peptide_parse_list)
      change_key(new_peptide['attr'], 'scan number', 'scan_id')
      change_key(new_peptide['attr'], 'retention time', 'retention_time')
      
      protein = protein_by_group_id[int(group_id)]
      i_source = i_sources[evidence['raw file']]
      protein['sources'][i_source]['peptides'].append(new_peptide)

  parse_proteins.count_peptides(proteins)
  
  return proteins, sources


def scale_log(r, max_v):
  v = math.log(r)
  scale_v = v/max_v
  if scale_v > 1.0:
    scale_v = 1.0
  if scale_v < -1.0:
    scale_v = -1.0
  return scale_v


def calculate_ratio_intensities(
    proteins, max_ratio=2.0, ratio_key='ratio h/l normalized'):
  for seqid in proteins.keys():
    for source in proteins[seqid]['sources']:
       for peptide in source['peptides']:
        ratio = peptide['attr'][ratio_key]
        if ratio is None or math.isnan(ratio) or ratio < 0:
          peptide['intensity'] = ""
        else:
          peptide['intensity'] = scale_log(ratio, max_ratio)


def calculate_lfq_ratio_intensities(
    proteins, experiment1, experiment2, max_ratio=2.0):
  for seqid in proteins:
    protein = proteins[seqid]
    peptide_intensities1 = []
    peptide_intensities2 = []
    for source in protein['sources']:
      peptides_by_seq = {}
      for peptide in source['peptides']:
        seq = peptide['sequence']
        peptides_by_seq.setdefault(seq, []).append(peptide)
      for seq, peptides in peptides_by_seq.items():
        intensities = {}
        for peptide in peptides:
          experiment = peptide['attr']['experiment']
          if not peptide['attr']['intensity']: 
            continue
          if experiment not in intensities:
            intensities[experiment] = []
          intensities[experiment].append(float(peptide['attr']['intensity']))
        if experiment1 in intensities and experiment2 in intensities:
          intensity1 = numpy.mean(intensities[experiment1])
          intensity2 = numpy.mean(intensities[experiment2])
          std1 = numpy.std(intensities[experiment1])
          std2 = numpy.std(intensities[experiment2])
          ratio = intensity1/intensity2
          std = ratio*numpy.sqrt((std1/intensity1)**2 + (std2/intensity2)**2)
          intensity = scale_log(ratio, max_ratio)
          std = numpy.round(std, 4)
          ratio = numpy.round(ratio, 4)
          peptide_intensities1.extend(intensities[experiment1])
          peptide_intensities2.extend(intensities[experiment2])
        elif experiment1 in intensities:
          ratio = float("inf")
          intensity = 2*max_ratio
          std = 0
          peptide_intensities1.extend(intensities[experiment1])
        elif experiment2 in intensities:
          ratio = 0.0
          intensity = -2*max_ratio
          std = 0
          peptide_intensities2.extend(intensities[experiment2])
        else:
          # neither experiment1 or experiment2 found in experiment column
          ratio = None
        for peptide in peptides:
          peptide['attr']['ratio'] = ratio
          peptide['intensity'] = intensity
          peptide['attr']['ratio_var'] = std
    sum2 = numpy.sum(peptide_intensities2)
    sum1 = numpy.sum(peptide_intensities1)
    if sum2 > 0.0:
      group_ratio = sum1/sum2
    else:
      group_ratio = float('inf')
    protein['attr']['ratio'] = group_ratio


if __name__ == '__main__':
  peptides, scans, protein_groups, evidence = read('../example/maxquant/silac')
  parse.save_data_dict(peptides, '../example/maxquant/peptides.dump')
  parse.save_data_dict(scans, '../example/maxquant/scans.dump')
  parse.save_data_dict(protein_groups, '../example/maxquant/protein_groups.dump')
  parse.save_data_dict(evidence, '../example/maxquant/evidence.dump')


