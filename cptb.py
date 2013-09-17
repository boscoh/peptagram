 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import os
import glob
from pprint import pprint

import peptagram.morpheus
import peptagram.mzml
import peptagram.proteins


def pepto_morpheus():
  proteins = peptagram.morpheus.get_proteins(
      'data/cptb/Morpheus/Results/MR20130807_C231PreSCX.protein_groups.tsv',
      'data/cptb/Morpheus/Results/MR20130807_C231PreSCX.PSMs.tsv',
      'data/cptb/Morpheus/modifications.tsv')
  peptagram.mzml.load_mzml(proteins, 0, 'data/cptb/MR20130807_C231PreSCX.mzML')
  peptagram.proteins.determine_unique_peptides(proteins)
  peptagram.proteins.count_peptides(proteins)
  out_dir = 'out/cptb'
  data = {
    'title': 'Corynbacterium Pseudo Tuberculosis', 
    'proteins': proteins,
    'source_labels': [''],
    'color_names': ['1.0', 'score/n', ''],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


if __name__ == '__main__':
  pepto_morpheus()
