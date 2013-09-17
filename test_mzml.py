 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint

import pyteomics.mzml


mzml = 'example/xtandem/Seq23282_E1O1.mzML'

indices = set()
ms_keys = set()
for spectrum in pyteomics.mzml.read(mzml):
  spectrum_id = spectrum['id']
  for key in spectrum.keys():
    if 'ms' in key.lower():
      print('{}: "{}"'.format(key, spectrum[key]))
  if spectrum_id in indices:
    print('Error same id', spectrum_id)
  else:
    indices.add(spectrum_id)
  if spectrum['ms level'] == 2:
    pprint(spectrum)
