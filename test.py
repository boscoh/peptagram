 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import math
import os
import json
import copy
import glob
import shutil
import logging


import peptagram.parse
import peptagram.maxquant
import peptagram.xtandem
import peptagram.morpheus
import peptagram.tpp
import peptagram.mzml
import peptagram.fasta
import peptagram.proteins



def clean_seqid(seqid):
  if '|' in seqid:
    return seqid.split('|')[1]
  else:
    return seqid


def test_xtandem():
  errors = [0.05, 0.025, 0.01]
  def get_protein():
    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'example/xtandem/interact.prot.xml', 
        'example/xtandem/interact.pep.xml', 
        n_peptide_cutoff=1, 
        is_skip_no_unique=True, errors=errors)
    for i in ['2', '3', '4']:
      fname = 'example/xtandem/Seq2328{}_E1O1.tandem'.format(i)
      for i, source_name in enumerate(source_names):
        basename = os.path.splitext(os.path.basename(fname))[0]
        if basename in source_name:
          i_source = i
          logging.debug('Matching pepxml source {} to xtandem {}'.format(basename, fname))
          break
      else:
        raise IOError('Couldn\'t match {} to {}'.format(basename, source_names))
      scans, fastas = peptagram.xtandem.read(fname)
      peptagram.proteins.load_fastas_into_proteins(proteins, fastas)
      peptagram.xtandem.load_scans_into_proteins(proteins, scans, i_source)
    return proteins
  proteins = peptagram.parse.memoize(
      get_protein, 'out/xtandem/proteins.dump', True)
  data = {
    'title': 'X!Tandem example',
    'proteins': proteins,
    'source_labels': ['yeast'],
    'color_names': ['P=1', 'P=0', ''],
    'mask_labels': map(str, errors),
  }
  peptagram.proteins.make_webapp_directory(data, 'out/xtandem')


def test_tpp():
  errors = [0.05, 0.025, 0.01]
  proteins, source_names = peptagram.tpp.get_proteins_and_sources(
      'example/tpp/hca-lysate-16.prot.xml', 
      'example/tpp/hca-lysate-16.pep.xml',
      errors=errors)
  peptagram.proteins.load_fasta_db_into_proteins(
      proteins, '../db/HUMAN.fasta')
  data = {
    'title': 'TPP example',
    'proteins': proteins,
    'source_labels': ['hca'],
    'color_names': ['P=1', 'P=0', ''],
    'mask_labels': map(str, errors),
  }
  peptagram.proteins.make_webapp_directory(data, 'out/tpp')


def test_maxquant():
  proteins, sources = peptagram.maxquant.get_proteins_and_sources(
      'example/maxquant/silac')
  peptagram.maxquant.calculate_ratio_intensities(proteins, max_ratio=1.5)
  peptagram.parse.save_data_dict(
      proteins, 
      'example/maxquant/silac/proteins.dump')
  peptagram.proteins.load_fasta_db_into_proteins(
      proteins, 
      'example/maxquant/silac/yeast_orf_trans_all_05-Jan-2010.fasta', 
      clean_seqid=clean_seqid)
  out_dir = 'out/maxq_silac/'
  data = {
    'title': 'Maxquant example', 
    'proteins': proteins,
    'source_labels': sources,
    'color_names': ['1.5', '1', '0.66'],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


def test_morpheus():
  proteins = peptagram.morpheus.get_proteins(
      'example/morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
      'example/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
      'example/morpheus//modifications.tsv'
      )
  peptagram.mzml.load_mzml(
      proteins, 0, 'example/morpheus/OK20130822_MPProtomap_KO1.mzML')
  peptagram.proteins.determine_unique_peptides(proteins)
  peptagram.proteins.count_peptides(proteins)
  peptagram.proteins.check_modifications(proteins)
  out_dir = 'out/morpheus'
  data = {
    'title': 'Morpheus Example', 
    'proteins': proteins,
    'source_labels': [''],
    'color_names': ['1.0', 'score/n', ''],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


if __name__ == '__main__':
  logging.basicConfig(level=logging.DEBUG)
  # test_tpp()
  # test_xtandem()
  # test_maxquant()
  test_morpheus()
  # todo: mascot and x!tandem get_proteins function
  