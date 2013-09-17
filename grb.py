 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import os
import glob
import copy
from pprint import pprint

import peptagram.morpheus
import peptagram.tpp
import peptagram.mzml
import peptagram.proteins


def flip_color_every_second_source(proteins):
  for protein in proteins.values():
    n_res = len(protein['sequence'])
    sources = protein['sources']
    n_source = len(sources)
    for i_source in range(0, n_source, 2):
      source1 = sources[i_source+1]
      for peptide in source1['peptides']:
        peptide['intensity'] *= -1


def morpheus():
  out_dir = 'out/grb/morpheus'
  protein_group_fname = 'data/Grb_Aug13/Morpheus/protein_groups.tsv'
  nums = range(1, 24)
  top_dir = 'data/Grb_Aug13/Morpheus/OK20130822_MPProtomap'
  modification_fname = 'data/Grb_Aug13/Morpheus/modifications.tsv'
  basenames = [top_dir + '1_wt', top_dir + '_KO']
  psm_fnames = ['{}{}.PSMs.tsv'.format(b, n) for n in nums for b in basenames]
  top_dir = '/Volumes/iago/Bosco/Grb_Aug13/mzML/OK20130822_MPProtomap'
  mzml_basenames = [top_dir + '1_wt', top_dir + '_KO']
  mzmls = ['{}{}'.format(b, n) for n in nums for b in mzml_basenames]
  mzmls = [m + '.mzML' for m in mzmls]
  sources = ['{}-{}'.format(e, n) for n in nums for e in ['WT', 'KO']]
  proteins = {}
  for psm_fname, mzml in zip(psm_fnames, mzmls):
    print("Reading", os.path.basename(psm_fname))
    these_proteins = peptagram.morpheus.get_proteins(
        protein_group_fname, psm_fname, modification_fname)
    print("Loading scans from ", os.path.basename(mzml))
    peptagram.mzml.load_mzml(these_proteins, 0, mzml)
    proteins = peptagram.proteins.merge_two_proteins(proteins, these_proteins)
  peptagram.proteins.determine_unique_peptides(proteins)
  peptagram.proteins.count_peptides(proteins)
  peptagram.proteins.check_modifications(proteins)
  flip_color_every_second_source(proteins)
  peptagram.parse.save_data_dict(proteins, out_dir + '/proteins.dump')
  data = {
    'title': 'Granzyme B', 
    'proteins': proteins,
    'source_labels': sources,
    'color_names': ['1.0', 'score/n', ''],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


def tpp():
  top_dir = 'data/Grb_Aug13/TPPXML/'
  fasta_db = top_dir + '../dbase/MOUSE_May13_clean.fasta'
  out_dir = 'out/grb/tpp'
  max_error = 0.01
  protxml = top_dir + 'CombinedALL/interactAll.prot.xml'
  protein_groups, protein_probs = peptagram.tpp.read_protxml(protxml)
  proteins = peptagram.tpp.make_proteins_from_protxml(protein_groups)
  basenames = [top_dir + 'interactGrB_wt', top_dir + 'interactGrB_KO']
  top_dir = '/Volumes/iago/Bosco/Grb_Aug13/mzML/OK20130822_MPProtomap'
  mzml_basenames = [top_dir + '1_wt', top_dir + '_KO']
  nums = range(1, 24)
  pepxmls = ['{}{}.pep.xml'.format(b, n) for n in nums for b in basenames]
  sources = ['{}-{}'.format(e, n) for n in nums for e in ['WT', 'KO']]
  mzmls = ['{}{}'.format(b, n) for n in nums for b in mzml_basenames]
  mzmls = [m + '.mzML' for m in mzmls]
  template_proteins = copy.deepcopy(proteins)
  for pepxml, mzml in zip(pepxmls, mzmls):
    print('Reading peptide-spectrum matches', os.path.basename(pepxml))
    dump_fname = os.path.basename(pepxml).replace('.pep.xml', '.dump')
    dump_fname = os.path.join(out_dir, dump_fname)
    scans_by_sources, peptide_probs = peptagram.parse.memoize(
        lambda: peptagram.tpp.read_pepxml(pepxml), dump_fname)
    these_proteins = copy.deepcopy(template_proteins)
    peptagram.tpp.load_pepxml(these_proteins, scans_by_sources)
    probability = 0.5
    print('probability cutoff', probability)
    peptagram.tpp.filter_peptides(these_proteins, probability)
    peptagram.mzml.load_mzml(these_proteins, 0, mzml)
    proteins = peptagram.proteins.merge_two_proteins(proteins, these_proteins)
  probability = peptagram.tpp.error_to_probability(protein_probs, max_error)
  peptagram.tpp.filter_proteins(proteins, probability)
  peptagram.proteins.load_fasta_db_into_proteins(proteins, fasta_db)
  peptagram.proteins.determine_unique_peptides(proteins)
  peptagram.proteins.count_peptides(proteins)
  peptagram.proteins.check_modifications(proteins)
  flip_color_every_second_source(proteins)
  peptagram.parse.save_data_dict(proteins, out_dir + '/proteins.dump')
  out_dir = 'out/grb/tpp'
  data = {
    'title': 'Granzyme B', 
    'proteins': proteins,
    'source_labels': sources,
    'color_names': ['1.0', 'prob', ''],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


def match_seqids(protein1, protein2):
  for seqid1 in protein1['attr']['seqids']:
    if seqid1 in protein2['attr']['seqids']:
      return True
  return False


def interweave(a_list):
  n = len(a_list)
  n_half = n/2
  result = []
  for i in range(n_half):
    result.append(a_list[i])
    result.append(a_list[i+n_half])
  return result


def compare_proteins():
  mor_proteins = eval(open('out/grb/morpheus/proteins.dump').read())
  tpp_proteins = eval(open('out/grb/tpp/proteins.dump').read())

  for proteins in [mor_proteins, tpp_proteins]:
    for seqid1, protein1 in proteins.items():
      protein1['attr']['seqid'] = seqid1
      protein1['attr']['seqids'] = [seqid1] + protein1['attr']['other_seqids']

  for seqid1, protein1 in tpp_proteins.items():
    sources1 = protein1['sources']
    n_source1 = len(sources1)
    print('---')
    print(protein1['attr']['seqids'], '->')
    match = False
    for i, (seqid2, protein2) in enumerate(mor_proteins.items()):
      if match_seqids(protein1, protein2):
        match = True
        break
    if match:
      sources2 = protein2['sources']
      n_source2 = len(sources2)
      print('  match:', protein2['attr']['seqids'])
      sources1.extend(sources2)
      if len(protein2['sequence']) > len(protein1['sequence']):
        protein1['sequence'] = protein2['sequence']
    else:
      sources2 = [{'peptides': []} for i in range(n_source2)]
      protein1['sources'] = sources1 + sources2
  for seqid2, protein2 in mor_proteins.items():
    for seqid1, protein1 in tpp_proteins.items():
      if match_seqids(protein1, protein2):
        break
    else:
      sources2 = mor_proteins[seqid2]['sources']
      sources1 = [{'peptides': []} for i in range(n_source1)]
      sources = sources1 + sources2
      tpp_proteins[seqid2] = protein2
      protein2['sources'] = sources
  
  for protein in proteins.values():
    n_res = len(protein['sequence'])
    sources = protein['sources']
    n_source = len(sources)
    protein_sequence = protein['sequence']
    for i_source in range(0, n_source):
      source = sources[i_source]
      peptides = source['peptides']
      for i_peptide in reversed(range(len(peptides))):
        peptide = peptides[i_peptide]
        peptide_sequence = peptide['sequence']
        i = protein_sequence.find(peptide_sequence)
        if i < 0:
          print("Error matching %s in %s" % (peptide_sequence, protein['attr']['seqid']))
          del peptides[i_peptide]
          continue
        peptide['i'] = i 
        peptide['j'] = i + len(peptide_sequence)
      for peptide in peptides:
        peptide['intensity'] = abs(peptide['intensity'])
        if i_source >= n_source/2:
          peptide['intensity'] = -peptide['intensity']
    protein['sources'] = interweave(protein['sources'])

  nums = range(1, 24)
  sources = ['{}{}-{}'.format(d, e, n) for n in nums for e in ['WT', 'KO'] for d in ['tpp', 'morpheus'] ]
  
  peptagram.proteins.determine_unique_peptides(proteins)
  peptagram.proteins.count_peptides(proteins)
  peptagram.proteins.check_modifications(proteins)

  out_dir = 'out/grb/compare'
  data = {
    'title': 'Granzyme B', 
    'proteins': tpp_proteins,
    'source_labels': sources,
    'color_names': ['1.0', 'prob', ''],
    'mask_labels': [],
  }
  peptagram.proteins.make_webapp_directory(data, out_dir)


def check_groups():
  mor_proteins = eval(open('out/grb/tpp/proteins.dump').read())
  groups = []
  for seqid, protein in mor_proteins.items():
    group = [seqid] + protein['attr']['other_seqids']
    groups.append(group)
  n_group = len(groups)
  print(n_group)
  for i_group1, group1 in enumerate(groups):
    intersect = []
    for i_group2 in range(i_group1+1, n_group):
      group2 = groups[i_group2]
      match = False
      for seqid in group1:
        if seqid in group2:
          match = True
          break
      if match:
        intersect.append(i_group2)
    if intersect:
      print(groups[i_group1])
      for i_group2 in intersect:
        print('  ', groups[i_group2])



if __name__ == '__main__':
  # morpheus()
  tpp()
  # compare_proteins()
  # check_groups()



