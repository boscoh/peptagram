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

import peptagram.parse
import peptagram.maxquant
import peptagram.morpheus
import peptagram.xtandem
import peptagram.mzml
import peptagram.tpp

this_dir = os.path.abspath(os.path.dirname(__file__))


def save_data_js(data, js_fname):
  f = open(js_fname, 'w')
  f.write('var data = \n')
  f.write(json.dumps(data, indent=None))
  f.close()


def read_fasta(fasta_db):
  """
  Returns a lsit of seqids encountered, and a dictionary
  of sequences for each seqid.
  """
  seqids = []
  seqid = None
  proteins = {}
  for line in open(fasta_db):
    if line.startswith(">"):
      words = line.split()
      seqid = words[0][1:]
      description = line[len(words[0]):].lstrip()
      seqids.append(seqid)
      proteins[seqid] = {
        'sequence': '',
        'description': description,
      }
      continue
    if seqid is not None:
      words = line.split()
      if words:
        proteins[seqid]['sequence'] += words[0]
  return seqids, proteins
  

def change_seqids_in_proteins(proteins, clean_seqid):
  seqids = proteins.keys()
  for seqid in seqids:
    new_seqid = clean_seqid(seqid)
    if seqid != new_seqid:
      proteins[new_seqid] = proteins[seqid]
      del proteins[seqid]
      if 'attr' in proteins[new_seqid]:
        if 'seqid' in proteins[new_seqid]['attr']:
          proteins[new_seqid]['attr']['seqid'] = new_seqid


def load_fasta_db(
    proteins, fasta_db, clean_seqid=None, iso_leu_isomerism=False):
  seqids, fastas = read_fasta(fasta_db)
  if clean_seqid:
    change_seqids_in_proteins(proteins, clean_seqid)
    change_seqids_in_proteins(fastas, clean_seqid)
  for seqid in proteins.keys():
    protein = proteins[seqid]
    if seqid not in fastas:
      print("Error: %s not found in fasta database" % seqid)
      del proteins[seqid]
      continue
    protein_sequence = fastas[seqid]['sequence']
    protein['description'] = fastas[seqid]['description']
    protein['sequence'] = protein_sequence
    if iso_leu_isomerism:
      protein_sequence = protein_sequence.replace("L", "I")
    n_peptide = 0
    for source in protein['sources']:
      peptides = source['peptides']
      for i_peptide in reversed(range(len(peptides))):
        peptide = peptides[i_peptide]
        peptide_sequence = peptide['sequence']
        if iso_leu_isomerism:
          peptide_sequence = peptide_sequence.replace("L", "I")
        i = protein_sequence.find(peptide_sequence)
        if i < 0:
          print("Error matching %s in %s" % (peptide_sequence, seqid))
          del peptides[i_peptide]
          continue
        peptide['i'] = i 
        peptide['j'] = i + len(peptide_sequence)
      n_peptide += len(peptides)
    if n_peptide == 0:
      # print("Deleted {}: no peptides found".format(seqid))
      del proteins[seqid]


def transfer_newer_files(in_dir, out_dir):
  print(in_dir, out_dir)
  for src in glob.glob(os.path.join(in_dir, '*')):
    dst = os.path.join(out_dir, os.path.basename(src))
    if not os.path.isfile(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
      shutil.copy(src, dst)


def make_webapp_directory(data, out_dir):
  # sanity checks
  proteins = data['proteins']
  for seqid, protein in proteins.items():
    for source in protein['sources']:
      peptides = source['peptides']
      for peptide in peptides:
        if 'mask' not in peptide:
          peptide['mask'] = 0
      peptides.sort(key=lambda peptide: len(peptide['sequence']))
      peptides.sort(key=lambda peptide: peptide['i'])
  if 'source_labels' not in data:
    data['source_labels'] = []
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
  save_data_js(data, os.path.join(out_dir, 'data.js'))
  transfer_newer_files(os.path.join(this_dir, 'templates/pepto'), out_dir)
  print("Made peptograph in ", end='')
  print(os.path.abspath(os.path.join(out_dir, 'index.html')))


def relpath(path):
  return os.path.relpath(path, os.getcwd())


def merge_proteins(proteins1, proteins2):
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
      sources1 = [{'peptides': []} for i in range(n_source1)]
      peptides = sources1 + sources2
      proteins1[seqid] = protein2
      proteins1[seqid]['sources'] = peptides
  for seqid in proteins1:
    peptides = proteins1[seqid]['sources']
    if len(peptides) == n_source1:
      peptides.extend([{'peptides': []} for i in range(n_source2)])
  return proteins1


def clean_seqid(seqid):
  if '|' in seqid:
    return seqid.split('|')[1]
  else:
    return seqid


def flip_color_every_second_source(proteins):
  for protein in proteins.values():
    n_res = len(protein['sequence'])
    sources = protein['sources']
    n_source = len(sources)
    for i_source in range(0, n_source, 2):
      source1 = sources[i_source+1]
      for peptide in source1['peptides']:
        peptide['intensity'] *= -1


def color_double_occupancy(proteins):
  def get_occupancy(source, n_res):
    occupancy = [False for i in range(n_res)]
    for peptide in source['peptides']:
      for i_res in range(peptide['i'], peptide['j']):
        occupancy[i_res] = True
    return occupancy

  def is_peptide_occupied(occupancy, peptide):
    for i_res in range(peptide['i'], peptide['j']):
      if occupancy[i_res] is False:
        return False
    return True

  for protein in proteins.values():
    n_res = len(protein['sequence'])
    sources = protein['sources']
    n_source = len(sources)
    for i_source in range(0, n_source, 2):
      source0 = sources[i_source]
      source1 = sources[i_source+1]
      occupancy0 = get_occupancy(source0, n_res)
      occupancy1 = get_occupancy(source1, n_res)
      double = [o0 and o1 for o0, o1 in zip(occupancy0, occupancy1)]
      for peptide in source0['peptides']:
        if is_peptide_occupied(double, peptide):
          peptide['intensity'] = 0.0
      for peptide in source1['peptides']:
        if is_peptide_occupied(double, peptide):
          peptide['intensity'] = 0.0


def get_numbered_ext_fnames(top_dir, ext, indices, prefix=''):
  save_dir = os.getcwd()
  os.chdir(top_dir)
  result = []
  for i in indices:
    tag = '{}*{:02}{}'.format(prefix, i, ext)
    fnames = glob.glob(tag)
    fnames.sort(key=lambda x: len(x))
    fnames = map(os.path.abspath, fnames)
    result.extend(fnames)
  os.chdir(save_dir)
  return result


def pepto_klk4_tpp_mzml():
  kl4_dir = '/Users/bosco/Projects/proteome/data/kl4/'
  kl4_dir = '/Volumes/iago/Bosco/KL4/'
  pepto_dir = os.path.join(kl4_dir, 'pepto')
  fasta_db = os.path.join(kl4_dir, 'HUMAN_May13.fasta')
  errors = [0.01]
  max_error = max(errors)
  experiment_dirs = [
    # 'LNCaP_Exp1',
    # 'LNCaP_Exp2',
    # 'LNCaP_Exp3',
    'WMPY1_Exp1',
    'WMPY1_Exp2',
    'WMPY1_Exp3',
    ]
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactLNCaP.prot.xml')
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactWMPY-1.prot.xml')

  for experiment_dir in experiment_dirs:
    tpp_dir = os.path.join(kl4_dir, 'TPP', experiment_dir)
    out_dir = os.path.join(pepto_dir, experiment_dir)
    if not os.path.isdir(out_dir):
      os.makedirs(out_dir)
    experiment_range = range(1, 20)

    print('Read Proteins', protxml)
    def get_proteins():
      protein_groups, protein_probs = peptagram.tpp.read_protxml(protxml)
      proteins = peptagram.tpp.make_proteins_from_protxml(protein_groups)
      return proteins, protein_probs
    proteins, protein_probs = peptagram.parse.memoize(
        get_proteins, pepto_dir + '/template_proteins.dump')

    # load pepxml PSM into proteins
    pepxmls = get_numbered_ext_fnames(tpp_dir, '.pep.xml', experiment_range)
    template_proteins = copy.deepcopy(proteins)
    for pepxml in pepxmls:
      print('Read Peptide-Spectrum Matches', os.path.basename(pepxml))
      dump_fname = os.path.basename(pepxml).replace('.pep.xml', '.dump')
      dump_fname = os.path.join(out_dir, dump_fname)
      def read_pepxml():
        scans_by_sources, peptide_probs = peptagram.tpp.read_pepxml(pepxml)
        return scans_by_sources, peptide_probs
      scans_by_sources, peptide_probs = peptagram.parse.memoize(read_pepxml, dump_fname)
      these_proteins = copy.deepcopy(template_proteins)
      peptagram.tpp.load_pepxml(these_proteins, scans_by_sources)
      probability = 0.5
      print('probability cutoff', probability)
      peptagram.tpp.filter_peptides(these_proteins, probability)
      probabilities = [peptagram.tpp.error_to_probability(peptide_probs, e) for e in errors]
      proteins = merge_proteins(proteins, these_proteins)

    load_fasta_db(proteins, fasta_db)

    peptagram.parse.determine_unique_peptides(proteins)

    # clean up empty proteins 
    probability = peptagram.tpp.error_to_probability(protein_probs, max_error)
    print('protein probability cutoff', probability)
    peptagram.tpp.filter_proteins(proteins, probability)
    peptagram.parse.count_peptides(proteins, is_skip_no_unique=True)
    flip_color_every_second_source(proteins)
    color_double_occupancy(proteins)

    source_labels = []
    for i in experiment_range:
      source_labels.extend([str(i), ''])

    data = {
      'title': 'KLK4 versus untreated', 
      'proteins': proteins,
      'source_labels': source_labels,
      'color_names': ['KLK4', '', 'CONTROL'],
      'mask_labels': map(str, errors),
    }
    make_webapp_directory(data, out_dir)


if __name__ == '__main__':
  # pepto_klk4_mq_lfq_mzml()
  pepto_klk4_tpp_mzml()
