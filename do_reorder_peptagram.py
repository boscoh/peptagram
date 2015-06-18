#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-


import os
import json
import webbrowser

import peptagram.proteins
from peptagram import parse

import Tkinter as tk
import tkFileDialog

import tkform


def load_data_jsonp(data_jsonp):
    with open(data_jsonp) as f:
      txt = f.read()
    i_bracket = txt.find('(')
    if i_bracket > 0:
        txt = txt[i_bracket+1:]
    i_bracket = txt.rfind(')')
    if i_bracket > 0:
        txt = txt[:i_bracket]
    return json.loads(txt)


class ResortPeptagramForm(tkform.Form):
  """
  Application window for Tkinter.
  """ 
  def __init__(self, width=700, height=800, parent=None):
    tkform.Form.__init__(self, parent, width, height)

    self.title('Peptagram Rearranger')

    self.push_text("Reorder Peptagram", 30)
    self.push_line()
    self.push_spacer()

    self.push_text("REQUIRED INPUT DATA", 16)
    self.push_text("Load existing peptagram directories; drag arrow to reorder; edit labels for peptagram")

    self.push_peptagram_loader(
        'in_peptagram', load_dir_text='+ peptagram directory')
    
    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("REORDERED-COMBINED PEPTAGRAM PARAMETERS", 16)
    self.push_labeled_param(
        'title', 'peptagram title', 'Reordered Peptagram')
    self.push_labeled_param(
        'out_dir', 'output directory', 'peptagram-reordered', load_dir_text='select')

    self.push_spacer()
    self.push_line()
    self.push_spacer()

    self.push_text("RESULTS", 16)
    self.push_submit()
    self.push_output()

  def make_source_label(self, i_data, i_source):
    return "Peptagram%d-Source%d" % (i_data+1, i_source+1)

  def extract_indices_from_source_label(self, source_label):
    s = source_label.replace('Peptagram', '').replace('Source', '')
    a, b = map(int, s.split('-'))
    return a-1, b-1

  def push_peptagram_loader(self, param_id, load_dir_text):
    self.file_list_loader = tkform.ReorderableList(self.interior)
    self.datas = []

    def load_peptagram():
      pep_dir = tkFileDialog.askdirectory(title=load_dir_text)
      try:
        self.print_output('Loading peptagram in ' + pep_dir + '...')
        data = load_data_jsonp(os.path.join(pep_dir, 'data.jsonp'))
      except:
        self.print_exception()

      self.datas.append(data)
      i_data = len(self.datas) - 1
      for i_source, label in enumerate(data['source_labels']):
        source_label = self.make_source_label(i_data, i_source)
        self.file_list_loader.add_entry_label(source_label, label)

    load_pep_button = tk.Button(self.interior, text=load_dir_text, command=load_peptagram)
    self.push_row(load_pep_button)

    self.push_row(self.file_list_loader)
    self.mouse_widgets.append(self.file_list_loader)
    self.param_entries[param_id] = self.file_list_loader

  def run(self, params):

    self.print_output("Building new peptagram...\n")
    seqids = []
    proteins = {}
    for data in self.datas:
      for seqid, protein in data['proteins'].iteritems():
        if seqid not in proteins:
          new_protein = peptagram.proteins.new_protein(seqid)
          new_protein['sources'] = []
          for key in ['sequence', 'description']:
            new_protein[key] = protein[key]
          new_protein['attr'] = protein['attr']
          seqids.append(seqid)
          proteins[seqid] = new_protein

    new_sources = []
    labels = []
    for source_label, label in params['in_peptagram']:
      i_data, i_source = self.extract_indices_from_source_label(source_label)
      new_sources.append((i_data, i_source))
      labels.append(label)

    for seqid in seqids:
      protein = proteins[seqid]
      for i_data, i_source in new_sources:
        old_data = self.datas[i_data]
        if seqid in old_data['proteins']:
          old_protein = old_data['proteins'][seqid]
          old_source = old_protein['sources'][i_source]
          protein['sources'].append(old_source)
        else:
          protein['sources'].append({'matches':[]})

    data = {
      'title': params['title'],
      'proteins': proteins,
      'source_labels': labels,
      'color_names': ['1.5', '1', '0.66'],
      'mask_labels': [],
    }

    out_dir = os.path.abspath(params['out_dir'])
    peptagram.proteins.make_graphical_comparison_visualisation(
        data, out_dir)
    self.print_output(
        'Successfully built peptagram webpage (%s):\n' % \
            parse.size_str(os.path.join(out_dir, 'data.js')))

    html = os.path.join(out_dir, 'index.html')
    cmd_fn = lambda: webbrowser.open('file://' + html)
    self.print_output(html, cmd_fn)


ResortPeptagramForm(800, -150).mainloop()




