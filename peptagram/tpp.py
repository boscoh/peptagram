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
    match['modified_sequence'] = match['peptide']

    match['other_seqids'] = []
    for alt_protein in search_hit_elem.findall(parse.fixtag('', 'alternative_protein', nsmap)):
      match['other_seqids'].append(alt_protein.attrib['protein'])

    match['modifications'] = []
    for modified_elem in search_hit_elem.findall(parse.fixtag('', 'modification_info', nsmap)):
      attr = parse.parse_attrib(modified_elem)
      match['modified_sequence'] = attr['modified_peptide']
      for modification_elem in modified_elem.findall(parse.fixtag('', 'mod_aminoacid_mass', nsmap)):
        attr = parse.parse_attrib(modification_elem)
        attr['i'] = attr['position'] - 1
        del attr['position']
        match['modifications'].append(attr)

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


class PepxmlReader(object):
  def __init__(self, pepxml):
    self.pepxml = pepxml
    self.distribution = None
    self.source_names = []
    self.i_source = None
    self.nsmap = {}
    self.is_debug = logger.root.level <= logging.DEBUG

  def __iter__(self):
    if self.is_debug:
      fname = self.pepxml + '.dump'
      logging.debug('Dumping pepxml reads into ' + fname)
      self.debug_file = open(fname, 'w')
      self.debug_file.write('[\n')
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
          for match in scan['matches']:
            fpe = probability_to_error(self.distribution, match['probability'])
            if fpe is None:
              print("WTF", match['probability'], self.distribution)
            else:
              match['fpe'] = fpe
          if self.i_source is not None:
            scan['source'] = self.source_names[self.i_source]
          if self.is_debug:
            pprint(scan, stream=self.debug_file)
            self.debug_file.write(',\n')
          yield scan
          elem.clear()
        elif elem.tag == parse.fixtag('', 'peptideprophet_summary', self.nsmap):
          self.distribution = parse_peptide_probabilities(elem, self.nsmap)
          if self.distribution[0]['prob'] < 1.0:
            self.distribution.insert(0, {'prob':1.0, 'error':0.0})
          if self.distribution[-1]['prob'] > 0.0:
            self.distribution.append({'prob':0.0, 'error':1.0})
          if self.is_debug:
            fname = self.pepxml + '.distribution.dump'
            pprint(self.distribution, open(fname, 'w'))
          elem.clear()
    if self.is_debug:
      self.debug_file.write(']\n')
      self.debug_file.close()


def make_peptide(pepxml_match, pepxml_scan, source):
  peptide = {
    'sequence': pepxml_match['peptide'],
    'modified_sequence': pepxml_match['modified_sequence'],
    'intensity': pepxml_match['probability'],
    'mask': pepxml_match['fpe'],
    'attr': {
      'pepxml_id': pepxml_scan['index'],
      'scan_id': pepxml_scan['start_scan'],
      'charge': pepxml_scan['assumed_charge'],
      'expect': pepxml_match['expect'],
      'modifications': pepxml_match['modifications'],
      'probability': pepxml_match['probability'],
      'missed_cleavages': pepxml_match['num_missed_cleavages'],
      'mass': pepxml_scan['precursor_neutral_mass'],
      'mass_diff': pepxml_match['massdiff'],
      'source': parse.basename(source),
    }
  }
  def grab_opt(peptide_key, scan_key, source_dict):
    if scan_key in source_dict:
      peptide['attr'][peptide_key] = source_dict[scan_key]
  grab_opt('retention_time', 'retention_time_sec', pepxml_scan)
  grab_opt('score', 'ionscore', pepxml_match)
  grab_opt('homology', 'homologyscore', pepxml_match)
  grab_opt('identity', 'identityscore', pepxml_match)
  peptide['attr']['matched_ions'] = str(pepxml_match['num_matched_ions'])
  peptide['attr']['matched_ions'] += '/'
  peptide['attr']['matched_ions'] += str(pepxml_match['tot_num_ions'])
  return peptide


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

    for parameter_elem in protein_elem.findall(parse.fixtag('', 'parameter', nsmap)):
      key = parameter_elem.attrib['name']
      val = parameter_elem.attrib['value']
      protein[key] = val

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
      peptide['modifications'] = []
      peptide['modified_sequence'] = peptide['peptide_sequence']
      for modified_elem in peptide_elem.findall(parse.fixtag('', 'modification_info', nsmap)):
        attr = parse.parse_attrib(modified_elem)
        peptide['modified_sequence'] = attr['modified_peptide']
        for modification_elem in modified_elem.findall(parse.fixtag('', 'mod_aminoacid_mass', nsmap)):
          attr = parse.parse_attrib(modification_elem)
          peptide['modifications'].append(attr)

    group['proteins'].append(protein)
  return group


class ProtxmlReader(object):
  def __init__(self, protxml):
    self.protxml = protxml
    self.distribution = None
    self.nsmap = {}
    self.debug_file = None
    self.is_debug = logger.root.level <= logging.DEBUG
    self.nsmap = {}

  def __iter__(self):
    if self.is_debug:
      fname = self.protxml + '.dump'
      logging.debug('Dumping protxml reads into ' + fname)
      self.debug_file = open(fname, 'w')
      self.debug_file.write('{\n')
    for event, elem in etree.iterparse(self.protxml, events=('end', 'start-ns')):
      if event == 'start-ns':
        self.nsmap.update({elem})
      if event == 'end':
        if elem.tag == parse.fixtag('', 'protein_group', self.nsmap):
          group = parse_protein_group(elem, self.nsmap)
          yield group
          if self.is_debug:
            pprint(group, stream=self.debug_file)
            self.debug_file.write(',\n')
          elem.clear()
        elif elem.tag == parse.fixtag('', 'proteinprophet_details', self.nsmap):
          self.distribution = parse_protein_probabilities(elem, self.nsmap)
          if self.is_debug:
            fname = self.protxml + '.distribution.dump'
            pprint(self.distribution, open(fname, 'w'))
          elem.clear()
    if self.is_debug:
      self.debug_file.write('}\n')
      self.debug_file.close()


def make_protein(protxml_protein):
  protein = {
    'description': protxml_protein['description'],
    'attr': { 
      'group_id': protxml_protein['group_number'],
      'length': protxml_protein['prot_length'],
      'seqid': protxml_protein['protein_name'],
      'sibling': protxml_protein['group_sibling_id'],
      'other_seqids': protxml_protein['other_seqids'],
      'probability': protxml_protein['probability'],
      'percent_coverage': '-',
    },
    'sources': [],
  }
  if 'percent_coverage' in protxml_protein:
    protein['attr']['percent_coverage'] = protxml_protein['percent_coverage']
  protein['protxml_peptides'] = protxml_protein['peptides']
  return protein


def generate_proteins_from_protxml(protxml):
  proteins = {}
  protxml_reader = ProtxmlReader(protxml)
  for protein_group in protxml_reader:
    for protxml_protein in protein_group['proteins']:
      seqid = protxml_protein['protein_name']
      if seqid in proteins:
        logger.warning("%s found in multiple protein groups" % seqid)
      proteins[seqid] = make_protein(protxml_protein)
  protein_probs = protxml_reader.distribution
  return proteins, protein_probs


def filter_proteins(proteins, prob_cutoff):
  for seqid in proteins.keys():
    protein = proteins[seqid]
    if protein['attr']['probability'] < prob_cutoff:
      del proteins[seqid]
      continue


def get_protein_by_seqid(proteins):
  protein_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    protein_by_seqid[seqid] = protein
    protein = proteins[seqid]
    for alt_seqid in protein['attr']['other_seqids']:
      protein_by_seqid[alt_seqid] = protein
  return protein_by_seqid


def add_source(proteins):
  for protein in proteins.values():
    protein['sources'].append({'peptides': []})


def load_pepxml(proteins, pepxml_fname, prob_cutoff=None, error_cutoff=None, source_names=[]):
  logging.debug('Peptide error cutoff: {}'.format(error_cutoff))
  protein_by_seqid = get_protein_by_seqid(proteins)
  n_source = len(proteins.values()[0]['sources'])
  pepxml_reader = PepxmlReader(pepxml_fname)
  source_name_indices = {}
  for scan in pepxml_reader:
    for match in scan['matches']:
      if error_cutoff is not None and match['fpe'] > error_cutoff:
        continue
      seqids = [match['protein']] + match['other_seqids']
      for seqid in seqids:
        if seqid not in protein_by_seqid:
          logger.warning('{} from scan {} not found in protxml'.format(seqid, scan['index']))
          continue
        protein = protein_by_seqid[seqid]
        if scan['source'] in source_name_indices:
          i_source = source_name_indices[scan['source']]
        else:
          add_source(proteins)
          n_source += 1
          i_source = n_source-1
          source_name_indices[scan['source']] = i_source
        peptides = protein['sources'][i_source]['peptides']
        scan_id = scan['start_scan']
        scan_ids = [p['attr']['scan_id'] for p in peptides]
        if scan_id not in scan_ids:
          peptide = make_peptide(match, scan, scan['source'])
          for protxml_peptide in protein['protxml_peptides']:
            if protxml_peptide['charge'] == scan['assumed_charge'] and \
                protxml_peptide['modified_sequence'] == match['modified_sequence']:  
              peptide['attr']['probability_protxml'] = protxml_peptide['nsp_adjusted_probability']
              peptide['attr']['is_contributing_evidence'] = protxml_peptide['is_contributing_evidence']
          peptides.append(peptide)
  source_names.extend(pepxml_reader.source_names)


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


def probability_to_error(distribution, prob):
  """
  Given a False-Positive-Error vs. Probability distribution,
  Cacluates the probability for an acceptable FPE.
  """
  fractionate = lambda a0, a1, a: (a-a0)/(a1-a0)
  interpolate = lambda a0, a1, f: a0 + f*(a1-a0)
  n = len(distribution)
  for i in range(1, n):
    prob0 = distribution[i-1]['prob']
    prob1 = distribution[i]['prob']
    if prob0 >= prob >= prob1:
      error0 = distribution[i-1]['error']
      error1 = distribution[i]['error']
      f = fractionate(prob0, prob1, prob)
      return interpolate(error0, error1, f)
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


def count_tpp_indep_spectra(proteins):
  for seqid in proteins:
    protein = proteins[seqid]
    n_unique_peptide = 0
    n_unique_spectrum = 0
    unique_sequence_set = set()
    for source in protein['sources']:
      for peptide in source['peptides']:
        sequence = str(peptide['attr']['charge']) + peptide['modified_sequence']
        if 'is_contributing_evidence' in peptide['attr'] and peptide['attr']['is_contributing_evidence'] == 'Y':
          n_unique_spectrum += 1
          unique_sequence_set.add(sequence)
    n_unique_peptide = len(unique_sequence_set)
    protein['attr']['n_indep_spectra'] = n_unique_spectrum 
  # calculate percentages of scans
  n_unique_spectrum_total = 0
  for seqid in proteins:
    protein = proteins[seqid]
    n_unique_spectrum_total += protein['attr']['n_indep_spectra']
  for seqid in proteins:
    protein = proteins[seqid]
    if n_unique_spectrum_total == 0:
      percent = 0
    else:
      percent = 100.0*protein['attr']['n_indep_spectra']/n_unique_spectrum_total
    protein['attr']['percent_indep_spectra'] = float('%.2f' % percent)


def get_proteins_and_sources(
    protxml, pepxmls, peptide_error=0.01, protein_error=0.01):
  """
  Returns a proteins dictionary and list of source names.
  """
  logger.info('Loading protxml ' + protxml)
  proteins, protein_probs = generate_proteins_from_protxml(protxml)

  source_names = []
  for pepxml in pepxmls:
    logger.info('Loading pepxml ' + pepxml)
    load_pepxml(proteins, pepxml, error_cutoff=peptide_error, source_names=source_names)
    
  count_tpp_indep_spectra(proteins)

  probability = error_to_probability(protein_probs, protein_error)
  filter_proteins(proteins, probability)

  if logger.root.level <= logging.DEBUG:
    dump = protxml.replace('prot.xml', 'proteins.dump')
    logger.debug('Dumping protein data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  return proteins, source_names


if __name__ == '__main__':
  logging.basicConfig(level=logging.DEBUG)
  proteins, sources = get_proteins_and_sources(
    '../example/tpp/hca-lysate-16.prot.xml', 
    ['../example/tpp/hca-lysate-16.pep.xml'])
  
