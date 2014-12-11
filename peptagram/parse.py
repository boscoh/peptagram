# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint
import re
import os
import json
import glob
import ntpath, posixpath, macpath


"""
Utility parsing functions that is useful for
different proteomics data parsers.
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
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return float(str(f)[:slen])
    
    
def basename(fname):
  for m in ntpath, posixpath, macpath:
    if m.sep in fname:
      return m.basename(fname)
  return fname


def parse_string(s):
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


# XML helper functions


def fixtag(ns, tag, nsmap):
  return '{' + nsmap[ns] + '}' + tag


def parse_attrib(elem):
  result = {}
  for key, value in elem.attrib.items():
    result[key] = parse_string(value)
  return result


def parse_name_value(elem):
  attrib = parse_attrib(elem)
  return { attrib['name']: attrib['value'] }


def save_data_dict(data, fname):
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
  size = 0
  for fname in fnames:
    if os.path.isdir(fname):
      sub_fnames = glob.glob(os.path.join(fname, '*'))
      size += size_of_fnames(*sub_fnames)
    else:
      size += os.path.getsize(fname)
  return size


def size_str(*fnames):
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


