# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import xml.etree.ElementTree as etree

import logging
logger = logging.getLogger('prophet')

import parse
import proteins as parse_proteins


"""
Parsers for ProteinProphet and PeptideProphet XML search results.

Main API entry:

  get_proteins_and_sources(
      protxml, 
      pepxmls, 
      peptide_error=0.01, 
      protein_error=0.01,
      good_expect=1E-8,
      cutoff_expect=1E-2)
  
  returns a dictionary that organizes peptide-spectrum-matches
  around proteins; and a list of sources.

"""


fractionate = lambda a0, a1, a: (a-a0)/(a1-a0)
interpolate = lambda a0, a1, f: a0 + f*(a1-a0)


def error_to_probability(distribution, error):
  """
  Given a False-Positive-Error vs. Probability distribution,
  Cacluates the probability for an acceptable FPE.
  """
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
      logger.debug('Dumping pepxml reads into ' + fname)
      self.debug_file = open(fname, 'w')
      self.debug_file.write('[\n')
    for event, elem in etree.iterparse(self.pepxml, events=('start', 'end', 'start-ns')):
      if event == 'start-ns':
        self.nsmap.update({elem})
      elif event == 'start':
        if elem.tag == self.search_tag('msms_run_summary'):
          fname = elem.attrib['base_name']
          self.source_names.append(fname)
          self.i_source = len(self.source_names) - 1
      elif event == 'end':
        if elem.tag == self.search_tag('spectrum_query'):
          scan = self.parse_scan(elem)
          for match in scan['matches']:
            fpe = probability_to_error(self.distribution, match['probability'])
            if fpe is None:
              logger.warning("WTF", match['probability'], self.distribution)
            else:
              match['fpe'] = fpe
          if self.i_source is not None:
            scan['source'] = self.source_names[self.i_source]
          if self.is_debug:
            pprint(scan, stream=self.debug_file)
            self.debug_file.write(',\n')
          yield scan
          elem.clear()
        elif elem.tag == self.search_tag('peptideprophet_summary'):
          self.parse_peptide_probabilities(elem)
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

  def search_tag(self, tag):
    return parse.fixtag('', tag, self.nsmap)

  def findall(self, elem, tag):
    return elem.findall(self.search_tag(tag))

  def find(self, elem, tag):
    return elem.find(self.search_tag(tag))

  def parse_peptide_probabilities(self, elem):
    # try with error_point
    error_points = self.findall(elem, 'error_point')
    if len(error_points) == 0:
      charge = 0
      for charge_elem in self.findall(elem, 'roc_error_data'):
        if charge_elem.attrib['charge'] == 'all':
          error_points = self.findall(charge_elem, 'error_point')
          break
    self.distribution = []
    for elem in error_points:
        attrib = parse.parse_attrib(elem)
        self.distribution.append({
          'error': attrib['error'],
          'prob': attrib['min_prob'],
        })
    self.distribution.sort(key=lambda d:d['error'])

  def parse_scan(self, scan_elem):
    scan = parse.parse_attrib(scan_elem)
    scan['matches'] = []
    for search_elem in self.findall(scan_elem, "search_result"):
      search_hit_elem = search_elem[0] 
      pepxml_match = parse.parse_attrib(search_hit_elem)
      pepxml_match['modified_sequence'] = pepxml_match['peptide']

      pepxml_match['other_seqids'] = []
      for alt_protein in self.findall(search_hit_elem, 'alternative_protein'):
        pepxml_match['other_seqids'].append(alt_protein.attrib['protein'])

      pepxml_match['modifications'] = []
      for modified_elem in self.findall(search_hit_elem, 'modification_info'):
        attr = parse.parse_attrib(modified_elem)
        pepxml_match['modified_sequence'] = attr['modified_peptide']
        for modification_elem in self.findall(modified_elem, 'mod_aminoacid_mass'):
          attr = parse.parse_attrib(modification_elem)
          attr['i'] = attr['position'] - 1
          del attr['position']
          pepxml_match['modifications'].append(attr)

      for score_elem in self.findall(search_hit_elem, 'search_score'):
        pepxml_match.update(parse.parse_name_value(score_elem))

      for analysis_elem in self.find(search_hit_elem, 'analysis_result'):
        if analysis_elem.tag == self.search_tag('peptideprophet_result'):
          pepxml_match.update(parse.parse_attrib(analysis_elem))
          for param_elem in analysis_elem[0]:
            pepxml_match.update(parse.parse_name_value(param_elem))

      scan['matches'].append(pepxml_match)

    return scan


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
      logger.debug('Dumping protxml reads into ' + fname)
      self.debug_file = open(fname, 'w')
      self.debug_file.write('{\n')
    for event, elem in etree.iterparse(self.protxml, events=('end', 'start-ns')):
      if event == 'start-ns':
        self.nsmap.update({elem})
      if event == 'end':
        if elem.tag == self.search_tag('protein_group'):
          group = self.parse_protein_group(elem)
          yield group
          if self.is_debug:
            pprint(group, stream=self.debug_file)
            self.debug_file.write(',\n')
          elem.clear()
        elif elem.tag == self.search_tag('proteinprophet_details'):
          self.parse_protein_probabilities(elem)
          if self.is_debug:
            fname = self.protxml + '.distribution.dump'
            pprint(self.distribution, open(fname, 'w'))
          elem.clear()
    if self.is_debug:
      self.debug_file.write('}\n')
      self.debug_file.close()

  def search_tag(self, tag):
    return parse.fixtag('', tag, self.nsmap)

  def findall(self, elem, tag):
    return elem.findall(self.search_tag(tag))

  def find(self, elem, tag):
    return elem.find(self.search_tag(tag))

  def parse_protein_probabilities(self, elem):
    self.distribution = []
    for data_point in self.findall(elem, 'protein_summary_data_filter'):
      attrib = parse.parse_attrib(data_point)
      self.distribution.append({
        'error': attrib['false_positive_error_rate'],
        'prob': attrib['min_probability'],
      })
    self.distribution.sort(key=lambda d:d['error'])

  def parse_protein_group(self, elem):
    group = parse.parse_attrib(elem)
    group['proteins'] = []
    for protein_elem in self.findall(elem, 'protein'):
      protein = parse.parse_attrib(protein_elem)
      protein['group_number'] = group['group_number']

      for parameter_elem in self.findall(protein_elem, 'parameter'):
        key = parameter_elem.attrib['name']
        val = parameter_elem.attrib['value']
        protein[key] = val

      annotation_elem = self.find(protein_elem, 'annotation')
      if annotation_elem is not None:
        protein['description'] = annotation_elem.attrib['protein_description']

      protein['other_seqids'] = []
      for alt_protein in self.findall(protein_elem, 'indistinguishable_protein'):
        protein['other_seqids'].append(alt_protein.attrib['protein_name'])

      protein['other_seqids'] = protein['other_seqids']
      protein['protein_name'] = protein['protein_name']

      protein['peptides'] = []
      n_unique_peptide = 0
      for peptide_elem in self.findall(protein_elem, 'peptide'):
        peptide = parse.parse_attrib(peptide_elem)
        protein['peptides'].append(peptide)
        peptide['modifications'] = []
        peptide['modified_sequence'] = peptide['peptide_sequence']
        for modified_elem in self.findall(peptide_elem, 'modification_info'):
          attr = parse.parse_attrib(modified_elem)
          peptide['modified_sequence'] = attr['modified_peptide']
          for modification_elem in self.findall(modified_elem, 'mod_aminoacid_mass'):
            attr = parse.parse_attrib(modification_elem)
            peptide['modifications'].append(attr)

      group['proteins'].append(protein)

    return group


def make_proteins_from_protxml(protxml):
  protxml_reader = ProtxmlReader(protxml)
  proteins = {}
  for protein_group in protxml_reader:
    for protxml_protein in protein_group['proteins']:
      seqid = protxml_protein['protein_name']
      if seqid in proteins:
        logger.warning("%s found in multiple protein groups" % seqid)
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
      proteins[seqid] = protein
  return proteins, protxml_reader.distribution


def get_protein_by_seqid(proteins):
  protein_by_seqid = {}
  for seqid in proteins:
    protein = proteins[seqid]
    protein_by_seqid[seqid] = protein
    protein = proteins[seqid]
    for alt_seqid in protein['attr']['other_seqids']:
      protein_by_seqid[alt_seqid] = protein
  return protein_by_seqid


def make_match(pepxml_match, pepxml_scan, source):
  match = {
    'sequence': pepxml_match['peptide'],
    'intensity': pepxml_match['probability'],
    'modifications': pepxml_match['modifications'],
    'mask': pepxml_match['fpe'],
    'attr': {
      'modified_sequence': pepxml_match['modified_sequence'],
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
  def grab_opt(match_key, scan_key, source_dict):
    if scan_key in source_dict:
      match['attr'][match_key] = source_dict[scan_key]
  grab_opt('retention_time', 'retention_time_sec', pepxml_scan)
  grab_opt('score', 'ionscore', pepxml_match)
  grab_opt('homology', 'homologyscore', pepxml_match)
  grab_opt('identity', 'identityscore', pepxml_match)
  match['attr']['matched_ions'] = str(pepxml_match['num_matched_ions'])
  match['attr']['matched_ions'] += '/'
  match['attr']['matched_ions'] += str(pepxml_match['tot_num_ions'])
  return match


def load_pepxml_into_proteins(
    proteins, 
    pepxml_fname, 
    prob_cutoff=None, 
    error_cutoff=None, 
    source_names=[],
    good_expect=1E-8,
    cutoff_expect=1E-2):

  logger.debug('Peptide error cutoff: {}'.format(error_cutoff))
  protein_by_seqid = get_protein_by_seqid(proteins)
  n_source = len(proteins.values()[0]['sources'])
  pepxml_reader = PepxmlReader(pepxml_fname)
  source_name_indices = {}
  for scan in pepxml_reader:
    for pepxml_match in scan['matches']:
      if error_cutoff is not None and pepxml_match['fpe'] > error_cutoff:
        continue
      seqids = [pepxml_match['protein']] + pepxml_match['other_seqids']
      for seqid in seqids:
        if seqid not in protein_by_seqid:
          logger.warning('{} from scan {} not found in protxml'.format(seqid, scan['index']))
          continue

        protein = protein_by_seqid[seqid]

        if scan['source'] in source_name_indices:
          i_source = source_name_indices[scan['source']]
        else:
          for protein in proteins.values():
            protein['sources'].append({'matches': []})
          n_source += 1
          i_source = n_source-1
          source_name_indices[scan['source']] = i_source

        matches = protein['sources'][i_source]['matches']
        scan_id = scan['start_scan']
        scan_ids = [m['attr']['scan_id'] for m in matches]

        if scan_id not in scan_ids:
          match = make_match(pepxml_match, scan, scan['source'])
          for protxml_peptide in protein['protxml_peptides']:
            if protxml_peptide['charge'] == scan['assumed_charge'] and \
                protxml_peptide['modified_sequence'] == pepxml_match['modified_sequence']:  
              match['attr']['probability_protxml'] = protxml_peptide['nsp_adjusted_probability']
              match['attr']['is_contributing_evidence'] = protxml_peptide['is_contributing_evidence']
          if pepxml_match['expect'] > cutoff_expect:
            continue
          match['intensity'] = \
              parse_proteins.calc_minus_log_intensity(
                  pepxml_match['expect'], good_expect, cutoff_expect)
          matches.append(match)

  source_names.extend(pepxml_reader.source_names)


def count_independent_spectra(proteins):
  for seqid in proteins:
    protein = proteins[seqid]
    n_unique_peptide = 0
    n_unique_spectrum = 0
    unique_sequence_set = set()
    for source in protein['sources']:
      for peptide in source['matches']:
        sequence = str(peptide['attr']['charge']) + peptide['attr']['modified_sequence']
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
    protxml, 
    pepxmls, 
    peptide_error=0.01, 
    protein_error=0.01,
    good_expect=1E-8,
    cutoff_expect=1E-2):
  """
  Returns a proteins dictionary and list of source names.
  """

  logger.info('Loading protxml ' + protxml)
  proteins, protein_probs = make_proteins_from_protxml(protxml)

  source_names = []
  for pepxml in pepxmls:
    logger.info('Loading pepxml ' + pepxml)
    load_pepxml_into_proteins(
      proteins, 
      pepxml, 
      error_cutoff=peptide_error, 
      source_names=source_names,
      good_expect=good_expect,
      cutoff_expect=cutoff_expect)
    
  count_independent_spectra(proteins)

  prob_cutoff = error_to_probability(protein_probs, protein_error)
  for seqid in proteins.keys():
    if proteins[seqid]['attr']['probability'] < prob_cutoff:
      del proteins[seqid]

  if logger.root.level <= logging.DEBUG:
    dump = protxml.replace('prot.xml', 'proteins.dump')
    logger.debug('Dumping protein data structure to ' + dump)
    parse.save_data_dict(proteins, dump)

  return proteins, source_names


  
