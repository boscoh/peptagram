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

import peptagram.parse
import peptagram.proteins
import peptagram.maxquant
import peptagram.morpheus
import peptagram.xtandem
import peptagram.mzml
import peptagram.tpp
import peptagram.fasta

this_dir = os.path.abspath(os.path.dirname(__file__))

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
      i_peptide = peptide['i']
      j_peptide = i_peptide + len(peptide['sequence'])
      for i_res in range(i_peptide, j_peptide):
        occupancy[i_res] = True
    return occupancy

  def is_peptide_occupied(occupancy, peptide):
    i_peptide = peptide['i']
    j_peptide = i_peptide + len(peptide['sequence'])
    for i_res in range(i_peptide, j_peptide):
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


def pepto_klk4_tpp_iago():
  kl4_dir = '/Users/bosco/Projects/proteome/data/kl4/'
  kl4_dir = '/Volumes/iago/Bosco/KL4/'
  pepto_dir = os.path.join(kl4_dir, 'pepto')
  fasta_db = os.path.join(kl4_dir, 'HUMAN_May13.fasta')
  errors = [0.01]
  max_error = max(errors)
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactLNCaP.prot.xml')
  experiment_dirs = [
    'LNCaP_Exp1',
    'LNCaP_Exp2',
    'LNCaP_Exp3',
    ]
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactWMPY-1.prot.xml')
  experiment_dirs = [
    'WMPY1_Exp1',
    'WMPY1_Exp2',
    'WMPY1_Exp3',
    ]
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactMutLNCap.prot.xml')
  experiment_dirs = [
    'mutLNCaP_Exp1',
    'mutLNCaP_Exp2',
    'mutLNCaP_Exp3',
    ]
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactMutWPMY-1All.prot.xml')
  experiment_dirs = [
    'mutWMPY1_Exp1',
    'mutWMPY1_Exp2',
    'mutWMPY1_Exp3',
    ]

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
      proteins = peptagram.proteins.merge_two_proteins(proteins, these_proteins)

    peptagram.proteins.load_fasta_db_into_proteins(proteins, fasta_db)

    peptagram.proteins.determine_unique_peptides(proteins)

    # clean up empty proteins 
    probability = peptagram.tpp.error_to_probability(protein_probs, max_error)
    print('protein probability cutoff', probability)
    peptagram.tpp.filter_proteins(proteins, probability)
    peptagram.proteins.count_peptides(proteins, is_skip_no_unique=True)
    flip_color_every_second_source(proteins)
    color_double_occupancy(proteins)

    source_labels = []
    for i in experiment_range:
      source_labels.extend([str(i), ''])

    data = {
      'title': 'KLK4 versus untreated', 
      'proteins': proteins,
      'source_labels': source_labels,
      'color_names': ['KLK4', '', 'mutant'],
      'mask_labels': map(str, errors),
    }
    peptagram.proteins.make_peptograph_directory(data, out_dir)


def pepto_klk4_tpp_local():
  kl4_dir = '/Users/bosco/Projects/peptagram/data/kl4/'
  pepto_dir = os.path.join(kl4_dir, 'pepto')
  fasta_db = os.path.join(kl4_dir, 'HUMAN_May13.fasta')
  errors = [0.01]
  max_error = max(errors)
  protxml = os.path.join(kl4_dir, 'TPP/protxml', 'interactLNCaP.prot.xml')
  experiment_dirs = [
    'WMPY1_Exp1',
    ]

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
      proteins = peptagram.proteins.merge_two_proteins(proteins, these_proteins)

    peptagram.proteins.load_fasta_db_into_proteins(proteins, fasta_db)

    peptagram.proteins.determine_unique_peptides(proteins)

    # clean up empty proteins 
    probability = peptagram.tpp.error_to_probability(protein_probs, max_error)
    print('protein probability cutoff', probability)
    peptagram.tpp.filter_proteins(proteins, probability)
    peptagram.proteins.count_peptides(proteins, is_skip_no_unique=True)
    flip_color_every_second_source(proteins)
    color_double_occupancy(proteins)

    source_labels = []
    for i in experiment_range:
      source_labels.extend([str(i), ''])

    data = {
      'title': 'KLK4 versus untreated', 
      'proteins': proteins,
      'source_labels': source_labels,
      'color_names': ['KLK4', '', 'mutant'],
      'mask_labels': map(str, errors),
    }
    peptagram.proteins.make_peptograph_directory(data, out_dir)


if __name__ == '__main__':
  logging.basicConfig(level=logging.INFO)
  # pepto_klk4_tpp_iago()
  pepto_klk4_tpp_local()
