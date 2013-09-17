 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import xml.etree.ElementTree as etree
import json
import logging
import os

import parse
import proteins as parse_proteins

"""
Parsers to read TPP output files: PEPXML and PROTXML.
"""

logger = logging.getLogger('tpp')


def parse_scan(scan_elem, nsmap):
  scan = parse.parse_attrib(scan_elem)
  scan['matches'] = []
  tag = lambda tag_id: parse.fixtag('', tag_id, nsmap)
  for search_elem in scan_elem.findall(parse.fixtag('', "search_result", nsmap)):
    search_hit_elem = search_elem[0] 
    match = parse.parse_attrib(search_hit_elem)
    match['modifications'] = []
    mod_tag = tag('modification_info')+'/'+tag('mod_aminoacid_mass')
    for mod_elem in search_hit_elem.findall(mod_tag):
      modification = parse.parse_attrib(mod_elem)
      modification['i'] = modification['position'] - 1
      del modification['position']
      match['modifications'].append(modification)
    for score_elem in search_hit_elem.findall(tag('search_score')):
      match.update(parse.parse_name_value(score_elem))
    for analysis_elem in search_hit_elem.find(parse.fixtag('', 'analysis_result', nsmap)):
      if analysis_elem.tag == parse.fixtag('', 'peptideprophet_result', nsmap):
        match.update(parse.parse_attrib(analysis_elem))
        for param_elem in analysis_elem[0]:
          match.update(parse.parse_name_value(param_elem))
    scan['matches'].append(match)
  return scan


def parse_peptide_probabilities(elem, nsmap):
  # try with error_point
  error_points = elem.findall(parse.fixtag('', 'error_point', nsmap))
  if len(error_points) == 0:
    charge = 0
    for charge_elem in elem.findall(parse.fixtag('', 'roc_error_data', nsmap)):
      if charge_elem.attrib['charge'] == 'all':
        error_points = charge_elem.findall(parse.fixtag('', 'error_point', nsmap))
        break
  probs = []
  for elem in error_points:
      attrib = parse.parse_attrib(elem)
      probs.append({
        'error': attrib['error'],
        'prob': attrib['min_prob'],
      })
  probs.sort(key=lambda d:d['error'])
  return probs


def read_pepxml(pepxml):
  nsmap = {}
  probs = []
  scan_sources = []
  for event, elem in etree.iterparse(pepxml, events=('start', 'end', 'start-ns')):
    if event == 'start-ns':
      nsmap.update({elem})
    elif event == 'start':
      if elem.tag == parse.fixtag('', 'msms_run_summary', nsmap):
        scan_source = {
          'scans': [],
          'filename': elem.attrib['base_name'],
        }
        scan_sources.append(scan_source)
    elif event == 'end':
      if elem.tag == parse.fixtag('', 'spectrum_query', nsmap):
        scan = parse_scan(elem, nsmap)
        scan_source['scans'].append(scan)
        elem.clear()
      elif elem.tag == parse.fixtag('', 'peptideprophet_summary', nsmap):
        probs = parse_peptide_probabilities(elem, nsmap)
        elem.clear()
  return scan_sources, probs


class PepxmlReader(object):
  def __init__(self, pepxml, prob_cutoff=None, error_cutoff=None):
    self.pepxml = pepxml
    self.probs = None
    self.source_names = []
    self.i_source = None
    self.prob_cutoff = None
    self.error_cutoff = None
    self.nsmap = {}

  def iter(self):
    for event, elem in etree.iterparse(self.pepxml, events=('start', 'end', 'start-ns')):
      if event == 'start-ns':
        self.nsmap.update({elem})
      elif event == 'start':
        if elem.tag == parse.fixtag('', 'msms_run_summary', self.nsmap):
          fname = elem.attrib['base_name']
          self.source_names.append(fname)
          self.i_source = len(self.source_names) - 1
      elif event == 'end':
        if elem.tag == parse.fixtag('', 'spectrum_query', self.nsmap):
          scan = parse_scan(elem, self.nsmap)
          if self.i_source is not None:
            scan['source'] = self.source_names[self.i_source]
          if self.prob_cutoff is None or scan['probability'] >= self.prob_cutoff:
            yield scan
          elem.clear()
        elif elem.tag == parse.fixtag('', 'peptideprophet_summary', self.nsmap):
          self.probs = parse_peptide_probabilities(elem, self.nsmap)
          if self.prob_cutoff is None and self.error_cutoff is not None:
            self.prob_cutoff = error_to_probability(self.probs, self.prob_cutoff)
          elem.clear()


def parse_protein_probabilities(elem, nsmap):
  probs = []
  for data_point in elem.findall(parse.fixtag('', 'protein_summary_data_filter', nsmap)):
    attrib = parse.parse_attrib(data_point)
    probs.append({
      'error': attrib['false_positive_error_rate'],
      'prob': attrib['min_probability'],
    })
  probs.sort(key=lambda d:d['error'])
  return probs


def parse_protein_group(elem, nsmap):
  group = parse.parse_attrib(elem)
  group['proteins'] = []
  for protein_elem in elem.findall(parse.fixtag('', 'protein', nsmap)):
    protein = parse.parse_attrib(protein_elem)
    protein['group_number'] = group['group_number']

    annotation_elem = protein_elem.find(parse.fixtag('', 'annotation', nsmap))
    if annotation_elem is not None:
      protein['description'] = annotation_elem.attrib['protein_description']

    protein['other_seqids'] = []
    for alt_protein in protein_elem.findall(parse.fixtag('', 'indistinguishable_protein', nsmap)):
      protein['other_seqids'].append(alt_protein.attrib['protein_name'])

    protein['other_seqids'] = protein['other_seqids']
    protein['protein_name'] = protein['protein_name']

    protein['peptides'] = []
    n_unique_peptide = 0
    for peptide_elem in protein_elem.findall(parse.fixtag('', 'peptide', nsmap)):
      peptide = parse.parse_attrib(peptide_elem)
      protein['peptides'].append(peptide)
      # if peptide['is_nondegenerate_evidence'] == 'Y':

    group['proteins'].append(protein)
  return group


def read_protxml(protxml):
  nsmap = {}
  distribution = {}
  protein_groups = {}
  for event, elem in etree.iterparse(protxml, events=('end', 'start-ns')):
    if event == 'start-ns':
      nsmap.update({elem})
    if event == 'end':
      if elem.tag == parse.fixtag('', 'protein_group', nsmap):
        group = parse_protein_group(elem, nsmap)
        protein_groups[group['group_number']] = group
        elem.clear()
      elif elem.tag == parse.fixtag('', 'proteinprophet_details', nsmap):
        distribution = parse_protein_probabilities(elem, nsmap)
        elem.clear()
  return protein_groups, distribution


def resort_matches_by_seq(scans_by_sources):
  for source in scans_by_sources:
    matches_by_seq = {}
    for scan in source['scans']:
      for peptide in scan['matches']:
        seq = peptide['peptide']
        match = { 'peptide': peptide, 'scan': scan }
        if seq not in matches_by_seq:
          matches_by_seq[seq] = []
        matches_by_seq[seq].append(match)
    source['matches_by_seq'] = matches_by_seq


def make_protein(protxml_protein):
  protein = {
    'description': protxml_protein['description'],
    'attr': { 
      'group_id': protxml_protein['group_number'],
      'seqid': protxml_protein['protein_name'],
      'sibling': protxml_protein['group_sibling_id'],
      'other_seqids': protxml_protein['other_seqids'],
      'probability': protxml_protein['probability'],
    }
  }
  protein['percent_coverage'] = None
  if 'percent_coverage' in protxml_protein:
    protein['percent_coverage'] = protxml_protein['percent_coverage']
  return protein


def make_peptide(pepxml_peptide, pepxml_scan, source):
  peptide = {
    'sequence': pepxml_peptide['peptide'],
    'attr': {
      'pepxml_id': pepxml_scan['index'],
      'scan_id': pepxml_scan['start_scan'],
      'expect': pepxml_peptide['expect'],
      'retention_time': pepxml_scan['retention_time_sec'],
      'modifications': pepxml_peptide['modifications'],
      'source': parse.basename(source),
    }
  }
  peptide['attr']['matched_ions'] = str(pepxml_peptide['num_matched_ions'])
  peptide['attr']['matched_ions'] += '/'
  peptide['attr']['matched_ions'] += str(pepxml_peptide['tot_num_ions'])
  peptide['attr']['probability'] = pepxml_peptide['probability']
  peptide['attr']['missed_cleavages'] = pepxml_peptide['num_missed_cleavages']
  peptide['attr']['mass'] = pepxml_scan['precursor_neutral_mass']
  peptide['attr']['mass_diff'] = pepxml_peptide['massdiff']
  peptide['intensity'] = pepxml_peptide['probability']
  return peptide


def not_in_peptides(new_peptide, peptides):
  for peptide in peptides:
    if peptide['sequence'] == new_peptide['sequence']:
      if peptide['attr']['scan_id'] == new_peptide['attr']['scan_id']:
        return False
  return True


def make_proteins_from_protxml(protein_groups):
  proteins = {}
  for protein_group in protein_groups.values():
    for protxml_protein in protein_group['proteins']:
      seqid = protxml_protein['protein_name']
      if seqid in proteins:
        logger.warning("%s found in multiple protein groups" % seqid)
      protein = make_protein(protxml_protein)
      proteins[seqid] = protein
      protein['sources'] = []
  return proteins


def filter_proteins(proteins, prob_cutoff):
  for seqid in proteins.keys():
    if proteins[seqid]['attr']['probability'] < prob_cutoff:
      del proteins[seqid]


def get_protein_by_seqid(proteins):
  protein_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    protein_by_seqid[seqid] = protein
    protein = proteins[seqid]
    for alt_seqid in protein['attr']['other_seqids']:
      protein_by_seqid[alt_seqid] = protein
  return protein_by_seqid


def load_pepxml(proteins, scans_by_sources):
  n_source = len(scans_by_sources)
  i_source_offset = None
  protein_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    protein_by_seqid[seqid] = protein
    protein = proteins[seqid]
    if i_source_offset is None:
      i_source_offset = len(protein['sources'])
    for i in range(n_source):
      protein['sources'].append({ 'peptides': [] })
    for alt_seqid in protein['attr']['other_seqids']:
      protein_by_seqid[alt_seqid] = protein
  for i_source, source in enumerate(scans_by_sources):
    for scan in source['scans']:
      for peptide in scan['matches']:
        seqid = peptide['protein']
        if seqid not in protein_by_seqid:
          logger.warning('{} from scan {} not found in protxml'.format(seqid, scan['index']))
          continue
        protein = protein_by_seqid[seqid]
        peptide = make_peptide(peptide, scan, source['filename'])
        i = i_source + i_source_offset
        peptides = protein['sources'][i]['peptides']
        peptides.append(peptide)


def error_to_probability(distribution, error):
  """
  Given a False-Positive-Error vs. Probability distribution,
  Cacluates the probability for an acceptable FPE.
  """
  fractionate = lambda a0, a1, a: (a-a0)/(a1-a0)
  interpolate = lambda a0, a1, f: a0 + f*(a1-a0)
  n = len(distribution)
  for i in range(1, n):
    error0 = distribution[i-1]['error']
    error1 = distribution[i]['error']
    if error0 <= error < error1:
      prob0 = distribution[i-1]['prob']
      prob1 = distribution[i]['prob']
      f = fractionate(error0, error1, error)
      return interpolate(prob0, prob1, f)
  return None


def filter_peptides(proteins, probability):
  seqids = proteins.keys()
  for seqid in seqids:
    n_peptide = 0
    for source in proteins[seqid]['sources']:
      peptides = source['peptides']
      for i_peptide in reversed(range(len(peptides))):
        if peptides[i_peptide]['attr']['probability'] < probability:
          del peptides[i_peptide]
        else:
          n_peptide += 1
    if n_peptide == 0:
      del proteins[seqid]


def make_mask(proteins, probabilities):
  "The probabilities are used to assign an 'i_mask' to each peptide."
  for i_mask, probability in enumerate(sorted(probabilities)):
    for seqid in proteins:
      for source in proteins[seqid]['sources']:
        for peptide in source['peptides']:
          if probability < peptide['attr']['probability']:
            peptide['mask'] = i_mask


def get_proteins_and_sources(
    protxml, pepxml, 
    n_peptide_cutoff=1, 
    is_skip_no_unique=True,
    errors = [0.01]):
  """
  Basic structure proteins in YAML formt.
    "sample_seqid": 
      sequence: "AAAAAAAAAA"
      description: "sample protein"
      attr:
        param: value
      sources:
        -
          peptides
            -
              sequence: "AAA"
              i: 0
              j: 3
              attr:
                is_unique: True
                param: value
  """
  max_error = max(errors)
  protein_groups, protein_probs = read_protxml(protxml)
  proteins = make_proteins_from_protxml(protein_groups)

  dump_dir = os.path.dirname(protxml)
  if logger.root.level <= logging.DEBUG:
    dump = os.path.join(dump_dir, 'protxml.dump')
    logger.debug('Dumping protxml data structure to ' + dump)
    parse.save_data_dict(protein_groups, dump)
    dump = os.path.join(dump_dir, 'proterror.dump')
    logger.debug('Dumping protein error distribution to ' + dump)
    parse.save_data_dict(protein_probs, dump)

  scans_by_sources, peptide_probs = read_pepxml(pepxml)

  if logger.root.level <= logging.DEBUG:
    dump = os.path.join(dump_dir, 'pepxml.dump')
    logger.debug('Dumping pepxml data structure to ' + dump)
    parse.save_data_dict(scans_by_sources, dump)
    dump = os.path.join(dump_dir, 'peperror.dump')
    logger.debug('Dumping peptide error distribution to ' + dump)
    parse.save_data_dict(peptide_probs, dump)

  source_names = [scans['filename'] for scans in scans_by_sources]
  load_pepxml(proteins, scans_by_sources)
  probability = error_to_probability(peptide_probs, max_error)
  filter_peptides(proteins, probability)
  probabilities = [error_to_probability(peptide_probs, e) for e in errors]
  make_mask(proteins, probabilities)
  probability = error_to_probability(protein_probs, max_error)
  filter_proteins(proteins, probability)
  parse_proteins.determine_unique_peptides(proteins)
  parse_proteins.count_peptides(proteins, n_peptide_cutoff, is_skip_no_unique)

  if logger.root.level <= logging.DEBUG:
    dump = os.path.join(dump_dir, 'proteins.dump')
    logger.debug('Dumping protein data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  return proteins, source_names


if __name__ == '__main__':
  logging.basicConfig(level=logging.DEBUG)
  proteins, sources = get_proteins_and_sources(
    '../example/tpp/hca-lysate-16.prot.xml', 
    '../example/tpp/hca-lysate-16.pep.xml')
  
