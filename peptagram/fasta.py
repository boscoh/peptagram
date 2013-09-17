# -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint


def read_fasta(fasta_db):
  """
  Returns a lsit of seqids encountered, and a dictionary
  of sequences for each seqid.
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
  
