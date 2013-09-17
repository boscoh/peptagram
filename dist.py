 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import math
import os
import json
import glob
import shutil
from pprint import pprint

import peptagram.parse
import peptagram.maxquant
import peptagram.morpheus
import peptagram.tpp
import peptagram.mzml

import plot
import pylab


errors = [0.05, 0.025, 0.01]
scans = eval(open('example/mascot/scans.dump').read())
scores = []
for scan in scans.values():
  for match in scan['matches']:
    scores.append(match['score'])

x_centers, bins = plot.count_bins(scores, 100, min(scores), max(scores))
plot.pylab_bars(x_centers, bins)
plot.pylab.savefig('pepdist.png')