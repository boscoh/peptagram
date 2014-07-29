#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import glob
import os
import sys
import traceback
import webbrowser
import os
import re

import Tkinter as tk
import tkFileDialog

import peptagram.xtandem
import peptagram.proteins


class VerticalScrolledFrame(tk.Frame):
  """
  Tkinter scrollable frame!
  - place widgets in the 'interior' attribute
  - construct and pack/place/grid normally
  - only allows vertical scrolling
  - adapted from http://stackoverflow.com/a/16198198
  """

  def __init__(self, parent, *args, **kw):
    tk.Frame.__init__(self, parent, *args, **kw)            

    # create a canvas object and a vertical scrollbar for scrolling it
    vscrollbar = tk.Scrollbar(self, orient=tk.VERTICAL)
    vscrollbar.pack(fill=tk.Y, side=tk.RIGHT, expand=tk.FALSE)
    self.canvas = tk.Canvas(
        self, bd=0, highlightthickness=0,
        yscrollcommand=vscrollbar.set)
    self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.TRUE)
    vscrollbar.config(command=self.canvas.yview)

    # reset the view
    self.canvas.xview_moveto(0)
    self.canvas.yview_moveto(0)

    # create a frame inside the canvas which will be scrolled with it
    self.interior = tk.Frame(self.canvas)
    self.interior_id = self.canvas.create_window(
        0, 0, window=self.interior, anchor=tk.NW)

    # track changes to canvas, frame and updates scrollbar
    self.interior.bind('<Configure>', self._configure_interior)
    self.canvas.bind('<Configure>', self._configure_canvas)
    self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

  def _on_mousewheel(self, event):
    self.canvas.yview_scroll(-1*(event.delta), "units")

  def _configure_interior(self, event):
    # update the scrollbars to match the size of the inner frame
    size = (self.interior.winfo_reqwidth(), self.interior.winfo_reqheight())
    self.canvas.config(scrollregion="0 0 %s %s" % size)
    if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
      # update the canvas's width to fit the inner frame
      self.canvas.config(width=self.interior.winfo_reqwidth())

  def _configure_canvas(self, event):
    if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
      # update the inner frame's width to fill the canvas
      self.canvas.itemconfigure(self.interior_id, width=self.canvas.winfo_width())


class FileEntry():
  def __init__(self, parent, fname, label):
    self.parent = parent
    self.fname = fname
    self.label_stringvar = tk.StringVar()
    self.label_stringvar.set(label)
    self.fname_widget = tk.Label(parent, text=self.fname)
    self.label_widget = tk.Entry(parent, textvariable=self.label_stringvar)
    self.delete_widget = tk.Label(parent, text="x")
    self.num_stringvar = tk.StringVar()
    self.num_stringvar.set('')
    self.num_widget = tk.Label(parent, textvariable=self.num_stringvar)

  def add_to_grid(self, j):
    self.num_stringvar.set(u'\u2195')
    self.j = j
    self.num_widget.grid(column=0,row=j,sticky='W')
    self.fname_widget.grid(column=1,row=j,sticky='W')
    self.label_widget = tk.Entry(self.parent, textvariable=self.label_stringvar)
    self.label_widget.grid(column=2,row=j,sticky='W')
    self.delete_widget.grid(column=3,row=j,sticky='W')
 
  def grid_forget(self):
    self.fname_widget.grid_forget()
    self.delete_widget.grid_forget()
    self.label_widget.destroy()
    self.num_widget.grid_forget()

  def in_y(self, event):
    y = event.y_root
    y0 = self.num_widget.winfo_rooty()
    y1 = self.num_widget.winfo_height() + y0
    return y0 <= y <= y1

  def in_grab_char(self, event):
    x, y = event.x_root, event.y_root
    y0 = self.num_widget.winfo_rooty()
    y1 = self.num_widget.winfo_height() + y0
    x0 = self.num_widget.winfo_rootx()
    x1 = self.num_widget.winfo_width() + x0
    return y0 <= y <= y1 and x0 <= x <= x1



class FileListLoader(tk.Frame):
  def __init__(self, parent):
    self.parent = parent
    tk.Frame.__init__(self, parent)
    self.grid()
    self.entries = []

  def add_fnames(self, fnames):
    for fname in fnames:
      label = os.path.basename(fname)
      xml_entry = FileEntry(self, fname, label)
      self.entries.append(xml_entry)
    self.clear_frame()
    self.build_frame()

  def clear_frame(self):
    for entry in self.entries:
      entry.grid_forget()

  def delete_entry(self, i):
    self.clear_frame()
    del self.entries[i]
    self.build_frame()

  def get_delete_callback(self, i):
    return lambda event: self.delete_entry(i)

  def build_frame(self):
    for i, entry in enumerate(self.entries):
      entry.add_to_grid(i)
      entry.delete_widget.bind(
          "<ButtonPress-1>", 
          self.get_delete_callback(i))

  def get_i_from_y(self, event):
    for i, entry in enumerate(self.entries):
      if entry.in_y(event):
        return i
    return -1

  def get_i_from_xy(self, event):
    for i, entry in enumerate(self.entries):
      if entry.in_grab_char(event):
        return i
    return -1

  def mouse_down(self, event):
    self.i_select = self.get_i_from_xy(event)
    if self.i_select == -1:
      return
    entry = self.entries[self.i_select]
    entry.num_widget.configure(background='#FF9999')
    entry.fname_widget.configure(background='#FF9999')

  def mouse_up(self, event):
    if self.i_select == -1:
      return
    entry = self.entries[self.i_select]
    entry.num_widget.configure(background='white')
    entry.fname_widget.configure(background='white')

  def mouse_drag(self, event):
    self.i_mouse_drag = self.get_i_from_y(event)
    if self.i_select == -1 or self.i_mouse_drag == -1:
      return
    if self.i_select == self.i_mouse_drag:
      return
    i, j = self.i_mouse_drag, self.i_select
    self.entries[i], self.entries[j] = self.entries[j], self.entries[i]
    self.clear_frame()
    self.build_frame()
    self.i_select = self.i_mouse_drag
    self.i_mouse_drag = -1


class HyperlinkManager:
  """
  http://effbot.org/zone/tkinter-text-hyperlink.htm
  """

  def __init__(self, text):
    self.text = text
    self.text.tag_config("hyper", foreground="blue", underline=1)
    self.text.tag_bind("hyper", "<Enter>", self._enter)
    self.text.tag_bind("hyper", "<Leave>", self._leave)
    self.text.tag_bind("hyper", "<Button-1>", self._click)
    self.reset()

  def reset(self):
    self.links = {}

  def add(self, action):
    # add an action to the manager.  returns tags to use in
    # associated text widget
    tag = "hyper-%d" % len(self.links)
    self.links[tag] = action
    return "hyper", tag

  def _enter(self, event):
    self.text.config(cursor="hand2")

  def _leave(self, event):
    self.text.config(cursor="")

  def _click(self, event):
    for tag in self.text.tag_names(tk.CURRENT):
      if tag[:6] == "hyper-":
        self.links[tag]()
        return


class LabeledEntry(tk.Frame):
  """
  Creates a frame that holds a label, a button and a text entry
  in a row.
  """
  def __init__(
      self, parent, text, entry_text='', 
      load_file_text=None,
      load_dir_text=None):
    self.parent = parent
    tk.Frame.__init__(self, parent)
    self.grid()

    self.stringvar = tk.StringVar()
    self.stringvar.set(entry_text)

    self.label = tk.Label(self, text=text)
    i_column = 0
    self.label.grid(column=i_column, row=0)

    i_column += 1
    self.button = None
    if load_file_text:
      self.button = tk.Button(
          self, text=load_file_text, command=self.load_file)
      self.button.grid(column=i_column, row=0)
      i_column += 1
    if load_dir_text:
      self.button = tk.Button(
          self, text=load_dir_text, command=self.load_dir)
      self.button.grid(column=i_column, row=0)
      i_column += 1

    self.entry = tk.Entry(self, textvariable=self.stringvar)
    self.entry.grid(column=i_column, row=0)

  def load_file(self):
    fname = tkFileDialog.askopenfilename()
    self.stringvar.set(fname)

  def load_dir(self):
    fname = tkFileDialog.askdirectory()
    self.stringvar.set(fname)


def size_str(*fnames):
  size = sum(map(os.path.getsize, fnames))
  if size < 1E6:
    return "%.3f MB" % (size/1E6)
  else:
    return "%.f MB" % (size/1E6)


def fix_list(tcl_list):
  if isinstance(tcl_list, list) or isinstance(tcl_list, tuple): 
    return tcl_list
  regex = r""" 
    {.*?}   # text found in brackets
    | \S+   # or any non-white-space characters 
  """
  tokens = re.findall(regex, tcl_list, re.X)
  # remove '{' from start and '}' from end of string
  return [re.sub("^{|}$", "", i) for i in tokens]


class PeptagramLoaderApp(tk.Tk):
  """
  Application window for Tkinter.
  """ 
  def __init__(self, parent):
    self.parent = parent
    tk.Tk.__init__(self, parent)
    self.geometry("700x800")

    self.vscroll_frame = VerticalScrolledFrame(self)
    self.vscroll_frame.pack(fill=tk.BOTH, expand=tk.TRUE)

    self.interior = self.vscroll_frame.interior
    self.interior.configure(bd=30)

    self.add_xml_button = tk.Button(
        self.interior, text='+ .tandem', command=self.load_xmls)

    self.file_list_loader = FileListLoader(self.interior)

    self.job_entry = LabeledEntry(
        self.interior, 'peptagram title', 'X!Tandem peptagram')
    self.outdir_entry = LabeledEntry(
        self.interior, 'output directory', 'peptagram', load_dir_text='select')
    self.n_peak_entry = LabeledEntry(
        self.interior, 'number of peaks for each PSM', '50')
    self.expectation_entry = LabeledEntry(
        self.interior, '(optional) expect cutoff for PSM (eg 0.01)')
    self.exclude_entry = LabeledEntry(
        self.interior, '(optional) filename of seqids to exclude',
        load_file_text='select')
    self.include_entry = LabeledEntry(
        self.interior, '(optional) filename of seqids to include',
        load_file_text='select')

    self.submit_button = tk.Button(
        self.interior,  text='submit', command=self.make_peptagram)
    self.output_text = tk.Text(self.interior, state=tk.DISABLED)

    self.i_row = 0
    self.render()

    self.bind('<Button-1>', self.file_list_loader.mouse_down) 
    self.bind('<B1-Motion>', self.file_list_loader.mouse_drag) 
    self.bind('<ButtonRelease-1>', self.file_list_loader.mouse_up) 
     
  def push_row(self, widget):
    self.i_row += 1
    widget.grid(row=self.i_row, column=0, sticky=tk.W)

  def render(self):
    def line():
      return tk.Canvas(self.interior, width=500, height=1, bg="#999999")

    def spacer(height=1):
      return tk.Label(self.interior, height=height)

    def heading(text, fontsize=12):
      return tk.Label(
          self.interior, font=('defaultFont', fontsize), text=text)

    self.push_row(heading("Create peptagram from X!Tandem", 30))
    self.push_row(line())
    self.push_row(spacer(height=1))

    self.push_row(heading("X!TANDEM FILES", 20))
    self.push_row(heading("drag arrow to reorder; edit labels for peptagram"))
    self.push_row(self.add_xml_button)
    self.push_row(self.file_list_loader)
    
    self.push_row(spacer(height=1))
    self.push_row(line())
    self.push_row(spacer(height=1))

    self.push_row(heading("PARAMETERS", 20))
    self.push_row(self.job_entry)
    self.push_row(self.outdir_entry)
    self.push_row(self.n_peak_entry)
    self.push_row(self.expectation_entry)
    self.push_row(self.exclude_entry)
    self.push_row(self.include_entry)

    self.push_row(spacer(height=1))
    self.push_row(line())
    self.push_row(spacer(height=1))

    self.push_row(heading("OUTPUT", 20))
    self.push_row(self.submit_button)
    self.push_row(self.output_text)
    
  def load_xmls(self):
    fnames = tkFileDialog.askopenfilenames(title="Open X!Tandem Files")
    # fix for Windows where askopenfilenames fails to format the list
    fnames = fix_list(fnames)
    self.file_list_loader.add_fnames(fnames)

  def get_params(self):
    params = {
      'tandems': [e.fname for e in self.file_list_loader.entries],
      'labels': [e.label_stringvar.get() for e in self.file_list_loader.entries],
      'out_dir': self.outdir_entry.stringvar.get(),
      'title': self.job_entry.stringvar.get(),
      'expect': None,
      'excluded_seqids': [],
      'include_seqids': [],
      'n_peak': 50,
    }
    if len(params['tandems']) == 0:
      raise ValueError('No .tandem files were selected.')
    
    expect_str = self.expectation_entry.stringvar.get()
    if expect_str:
      params['expect'] = float(expect_str)

    excluded_fname = self.exclude_entry.stringvar.get()
    if excluded_fname:
      params['excluded_seqids'] = open(excluded_fname).read().split()

    included_fname = self.include_entry.stringvar.get()
    if included_fname:
      params['included_seqids'] = open(included_fname).read().split()

    n_peak = self.n_peak_entry.stringvar.get()
    params['n_peak'] = int(n_peak)
    
    return params

  def generate_peptagram_from_xtandem(self, params):
    proteins = {}
    for tandem in params['tandems']:
      s = "Processing " + tandem + " (%s)...\n" % size_str(tandem)
      self.output_text.insert(tk.INSERT, s)
      self.update()
      these_proteins = peptagram.xtandem.create_proteins_from_xtandem(
          tandem, 
          n_peak=params['n_peak'],
          expect_upper_cutoff=params['expect'],
          excluded_seqids=params['excluded_seqids'])
      proteins = peptagram.proteins.merge_two_proteins(
          proteins, these_proteins)
    data = {
      'title': params['title'],
      'proteins': proteins,
      'source_labels': params['labels'],
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': [],
    }
    peptagram.proteins.make_graphical_comparison_visualisation(
        data, params['out_dir'])

  def make_peptagram(self):
    self.output_text.configure(state=tk.NORMAL)
    self.output_text.delete(1.0, tk.END)
    self.update()
    try:
      params = self.get_params()
      self.output_text.insert(
          tk.INSERT, 
          'Generating peptagram from %s of files...\n' %
              size_str(*params['tandems']))
      self.update()
      self.generate_peptagram_from_xtandem(params)
      out_dir = os.path.abspath(params['out_dir'])
      self.output_text.delete(1.0, tk.END)
      self.output_text.insert(
          tk.INSERT, 
          'Successfully built peptagram webpage (%s):\n' % \
              size_str(os.path.join(out_dir, 'data.js')))
      link = HyperlinkManager(self.output_text)
      html = os.path.join(out_dir, 'index.html')
      callback = link.add(lambda: webbrowser.open('file://' + html))
      self.output_text.insert(tk.INSERT, html, callback)
    except:
      s = "\nTHERE WERE ERROR(S) IN PROCESSING THE PYTHON.\n"
      s += "Specific error described in the last line:\n\n"
      s += traceback.format_exc()
      self.output_text.insert(tk.INSERT, s)
    self.output_text.configure(state=tk.DISABLED)


if __name__ == "__main__":
    app = PeptagramLoaderApp(None)
    app.title('my application')
    app.mainloop()




