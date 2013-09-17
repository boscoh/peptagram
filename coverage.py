import os
import pprint
import json
import shutil
import glob

import envoy 

import peptagram.tpp
import peptagram.xtandem
import peptagram.morpheus
import peptagram.tpp_xtandem
import peptagram.parse

from jinja2 import Environment, FileSystemLoader
from hamlpy.ext import HamlPyExtension

template_dir = os.path.join(
    os.path.dirname(__file__), 'templates', 'coverage')



def sequence_coverage_html(sequence, peptides):
  def get_peptide(i):
    for peptide in peptides:
      if peptide['i'] <= i < peptide['j']:
        return peptide
    return None
  s = ""
  current_peptide = None
  for i in range(len(sequence)):
    peptide = get_peptide(i)
    if peptide:
      a_tag = "<a class='highlight' href='{}'>"
      a_tag = a_tag.format(peptide['link'])
    if peptide and current_peptide is None:
      s += a_tag
      current_peptide = peptide
    if not peptide and current_peptide:
      s += "</a>"
      current_peptide = None
    if i % 50 == 0:
      if current_peptide:
        s += "</a>"
      if i > 0:
        s += "<br/>"
      s += str(i).rjust(5) + ' '
      if current_peptide:
        s += a_tag
    elif i % 10 == 0:
      s += " "
    s += sequence[i]
  if current_peptide:
    s += "</a>"
  return s


def make_html_from_haml(haml, env, out_dir):
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
  html = out_dir + '/index.html'
  dirname, basename = os.path.split(haml)
  jinja2_env = Environment(
    loader=FileSystemLoader(dirname), 
    extensions=[HamlPyExtension])
  template =  jinja2_env.get_template(basename)
  html_text = template.render(env)
  with open(html, 'w') as f:
    f.write(html_text)
  print("Made", html)
  for f in glob.glob(template_dir + '/*js'):
    shutil.copy(f, out_dir)


def convert_proteins_to_html(proteins, out_dir, attr_keys):
  entries = []
  data = {}
  seqids = proteins.keys()
  for i_seqid, seqid in enumerate(sorted(seqids)):
    protein = proteins[seqid]
    sequence = protein['sequence']
    peptides = []
    for source in protein['sources']:
      for source_peptide in source['peptides']:
        peptide = {}
        for key in attr_keys:
          if key in source_peptide['attr']:
            peptide[key] = source_peptide['attr'][key]
          else:
            peptide[key] = ''
        peptide['sequence'] = source_peptide['sequence']
        peptide['i'] = source_peptide['i']
        peptide['j'] = source_peptide['j']
        if 'spectrum' in source_peptide:
          peptide['spectrum'] = source_peptide['spectrum']
        else:
          peptide['spectrum'] = []
        peptides.append(peptide)

      peptides.sort(key=lambda peptide: peptide['i'])    
      for i, peptide in enumerate(peptides):
        peptide['#'] = i+1
        peptide['id'] = seqid + '-' + str(i+1)
        peptide['link'] = '#' + peptide['id']
        data[peptide['link']] = peptide['spectrum']

      sequence_html = sequence_coverage_html(sequence, peptides)

      protein = {
        'description': protein['description'],
        'i_seqid': i_seqid + 1,
        'seqid': seqid,
        'seqid_link': '#' + seqid,
        'peptides': peptides,
        'sequence_html': sequence_html,
      }
      entries.append(protein)

  env = { 
    'title': 'This is an Example',
    'entries': entries, 
    'attr_keys': attr_keys,
    'data': data
  }

  haml = os.path.join(template_dir, 'template.haml')
  make_html_from_haml(haml, env, out_dir)



def test_showit_xtandem():
  html = 'peptide.html'
  protxml = 'example/xtandem/interact.prot.xml'
  pepxml = 'example/xtandem/interact.pep.xml'
  xtandem_xml = 'example/xtandem/Seq23282_E1O1.tandem'

  is_skip_no_unique = True
  is_only_one_sibling = True
  n_peptide_cutoff = 1
  attr_keys = [
    '#',
    'scan_id',
    'i',
    'sequence',
    'nsp_probability',
    'probability',
    'is_unique',
    'matched_ions',
    'missed_cleavages',
  ]
  dump = 'example/xtandem/proteins.dump'
  proteins = peptagram.tpp_xtandem.get_proteins(
      protxml, pepxml, xtandem_xml, errors=[0.05])
  open(dump, 'w').write(repr(proteins))
  proteins = eval(open(dump).read())
  if os.path.isfile(html): os.remove(html)
  convert_proteins_to_html(proteins, 'out/coverage_x/', attr_keys)


def test_showit_morpheus():
  morpheus = 'example/morpheus0/OK20130628_PlGQuick'
  proteins = peptagram.morpheus.read_proteins(morpheus)
  attr_keys = [
    'i',
    'sequence',
    'morpheus score',
    'matching products',
    'total products',
    'tot_num_ions',
    'scan number',
    'missed cleavages',
  ]
  convert_proteins_to_html(proteins, 'out/coverage_m', attr_keys)


if __name__ == "__main__":
  test_showit_morpheus()
  # test_showit_xtandem()
 
