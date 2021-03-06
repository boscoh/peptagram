# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint

import re
import os
import json
import glob
import logging
import ntpath, posixpath, macpath


"""
Utility parsing functions for strings and files.
"""


# string parsers

float_regex_pattern = r"""
^ 
[-+]? # optional sign
(?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
)
# followed by optional exponent part if desired
(?: [Ee] [+-]? \d+ ) ?
$
"""
float_regex = re.compile(float_regex_pattern, re.VERBOSE)


def round_decimal(f, n):
  "Truncates/pads a float f to n decimal places without rounding"
  slen = len('%.*f' % (n, f))
  return float(str(f)[:slen])
    
    
def basename(fname):
  "Returns the fname with the directory part removed"
  for m in ntpath, posixpath, macpath:
    if m.sep in fname:
      return m.basename(fname)
  return fname


def parse_string(s):
  "Converts a string to a float or int if matches numerical pattern"
  if re.search(r'^[-+]?\d+$', s):
    return int(s)
  elif float_regex.match(s):
    return float(s)
  else:
    return s


# tsv helper function

def splitter(s, convert_fn=None, delimiter=";"):
  s = str(s).strip()
  if not s:
    return []
  words = s.split(delimiter)
  if convert_fn:
    words = map(convert_fn, words)
  return words


def split_tab(line):
  if line[-1] == '\n':
    line = line[:-1]
  return line.split('\t')


def read_tsv(tsv_txt):
  """
  Reads a title top TSV file, converts all keys to lower case
  and returns a list of dictionaries.
  """
  f = open(tsv_txt, "UR")
  titles = [w.lower() for w in split_tab(f.readline())]
  results = []
  for line in f.readlines():
    group = {}
    for key, val in zip(titles, split_tab(line)):
      group[key] = parse_string(val)
    yield group
    del group


# XML helper functions


def fixtag(ns, tag, nsmap):
  return '{' + nsmap[ns] + '}' + tag


def parse_attrib(elem):
  "Returns a Python dictionary from an xml attrib dict"
  result = {}
  for key, value in elem.attrib.items():
    result[key] = parse_string(value)
  return result


def parse_name_value(elem):
  attrib = parse_attrib(elem)
  return { attrib['name']: attrib['value'] }


def save_data_dict(data, fname):
  "Prints out a Python dictionary to a formatted text file"
  with open(fname, 'w') as f:
    pprint(data, stream=f, indent=2)


def cache_result_to_json(fname):
  def actual_decorator(fn):
    def wrapped_fn(*args, **kwargs):
      if not os.path.isfile(fname):
        result = fn(*args, **kwargs)
        json.dump(result, open(fname, 'w'))
      else:
        result = json.load(open(fname))
      return result
    return wrapped_fn
  return actual_decorator


def memoize(fn, fname, force=False):
  if not force and os.path.isfile(fname):
    return eval(open(fname).read())
  else:
    result = fn()
    dirname = os.path.dirname(fname)
    if dirname is not '' and not os.path.isdir(dirname):
      os.makedirs(dirname)
    pprint(result, stream=open(fname, 'w'), width=80)
    return result


def size_of_fnames(*fnames):
  "Returns the size of the fnames that are the arguments of the fn"
  size = 0
  for fname in fnames:
    if os.path.isdir(fname):
      sub_fnames = glob.glob(os.path.join(fname, '*'))
      size += size_of_fnames(*sub_fnames)
    elif os.path.isfile(fname):
      size += os.path.getsize(fname)
  return size


def check_fnames(*fnames):
  for fname in fnames:
    fname = os.path.abspath(fname)
    if os.path.isfile(fname):
      continue
    elif os.path.isdir(fname):
      continue
    else:
      raise IOError(
          "%s does not exist" % fname)


def size_str(*fnames):
  "Formated size string that is more human readable"
  check_fnames(*fnames)
  size = size_of_fnames(*fnames)
  if size < 1E6:
    return "%.3f MB" % (size/1E6)
  else:
    return "%.f MB" % (size/1E6)


def read_word_file(fname):
  if not fname:
    return []
  if not os.path.isfile(fname):
    raise IOError("Couldn't find file: " + fname)
  text = open(fname).read()
  return text.split()


def read_fasta(fasta_db):
  """
  Returns list of seqids and a dictionary of
  proteins sequences as an attribute dictionary.
  """
  seqids = []
  seqid = None
  proteins = {}
  for line in open(fasta_db):
    if line.startswith(">"):
      words = line.split()
      seqid = words[0][1:]
      description = line[len(words[0]):].lstrip()
      seqids.append(seqid)
      proteins[seqid] = {
        'sequence': '',
        'description': description,
      }
      continue
    if seqid is not None:
      words = line.split()
      if words:
        proteins[seqid]['sequence'] += words[0]
  return seqids, proteins
  

class DictListWriter():
  """
  Writes a list of pprint-formatted dictionaries for debugging.
  Can be switched on/off by the is_debug flag.
  """
  def __init__(self, is_debug, fname):
    self.fname = fname
    self.is_debug = is_debug
    if self.is_debug:
      logging.debug('Dumping dict list to ' + self.fname)
      self.file = open(self.fname, 'w')
      self.file.write('[\n')

  def dump_dict(self, data_dict):
    if self.is_debug:
      pprint(data_dict, stream=self.file)
      self.file.write(',\n')

  def close(self):
    if self.is_debug:
      self.file.write(']\n')
      self.file.close()




