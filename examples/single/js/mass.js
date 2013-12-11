

atom_mass = {
  'H': 1.007825032,
  'O': 15.994914622,
}

ave_atom_mass = {
  'H': 1.00794,
  'O': 15.9994,
}

aa_monoisotopic_mass = {
  'A': 71.037114,
  'B': 110.000000,
  'C': 103.009184,
  'D': 115.026943,
  'E': 129.042593,
  'F': 147.068414,
  'G': 57.021464,
  'H': 137.058912,
  'I': 113.084064,
  'J': 110.000000,
  'K': 128.094963,
  'L': 113.084064,
  'M': 131.040485,
  'N': 114.042927,
  'O': 110.000000,
  'P': 97.052764,
  'Q': 128.058578,
  'R': 156.101111,
  'S': 87.032028,
  'T': 101.047678,
  'U': 110.000000,
  'V': 99.068414,
  'W': 186.079313,
  'X': 110.000000,
  'Y': 163.063329,
  'Z': 110.000000
}

ion_type_props = {
  'b': { 
    'atoms_to_add': [], // Cterm O+ loses electron to form triple bond to C
    'atoms_to_sub': [],
    'is_nterm': true, 
  },
  'y': { 
    'atoms_to_add': ['H', 'H', 'O'], // H+ on Nterm NH, OH on Cterm CO
    'atoms_to_sub': [],
    'is_nterm': false,
  },
}

/* returns a list of masses with len(sequence)+2.
  First position is for N-termius modifications, and
  last position is for C-terminus modifications */
function parse_terminal_aa_masses(sequence, aa_mass, modified_masses) {
  masses = [0.0];
  for (var i=0; i<sequence.length; i++) {
    var aa = sequence[i];
    var mass = aa_mass[aa];
    masses.push(mass);
  }
  masses.push(0.0);
  for (var j=0; j<modified_masses.length; j++) {
    var modified = modified_masses[j];
    var i_with_terminii = modified.i+1;
    if (i_with_terminii < masses.length) {
      masses[i_with_terminii] = Number(modified.mass);
    } else {
      console.log('warning modification.i longer than sequence', sequence);
    }
  }
  return masses;
}


function ion_label(i, charge, ion_type) {
  if (charge == 1) {
    return ion_type + i;
  } else {
    return ion_type + i + '(' + charge + '+)';
  }
}


function calculate_mz(aa_masses, charge, ion_type) {
  var mass = 0;
  for (var i=0; i<aa_masses.length; i++) {
    mass += aa_masses[i];
  }
  var atoms = ion_type_props[ion_type].atoms_to_add;
  for (var i=0; i<atoms.length; i++) {
    mass += atom_mass[atoms[i]]
  }  
  var atoms = ion_type_props[ion_type].atoms_to_sub;
  for (var i=0; i<atoms.length; i++) {
    mass -= atom_mass[atoms[i]]
  }  
  mass += atom_mass['H']*charge;
  if (charge == 1) {
    return mass;
  } 
  return mass/charge;
}


function calculate_peaks(terminii_aa_masses, charge, ion_type) {
  peaks = [];
  var n_all = terminii_aa_masses.length; 
  var n_seq = n_all - 2;
  for (var i_cut=0; i_cut<n_seq; i_cut++) {
    if (ion_type_props[ion_type].is_nterm) {
      fragment = terminii_aa_masses.slice(0,i_cut+2);
    } else {
      fragment = terminii_aa_masses.slice(i_cut+1,n_all);
    }
    var fragment_length = fragment.length-1
    if (fragment_length > 0) {
      mz = calculate_mz(fragment, charge, ion_type);
      // round to 4 decimal places
      mz = Math.round(1000*mz)/1000
      label = ion_label(fragment_length, charge, ion_type);
      peaks.push([mz, label]);
    }
  }
  return peaks;
}


function map_matched_ions(ion_type, sequence, peaks, mz_delta, modified_aa_masses, aa_mass) {
  var terminii_aa_masses = parse_terminal_aa_masses(sequence, aa_mass, modified_aa_masses);
  if (ion_type[0] == 'b') {
    var ion = 'b';
  } else {
    var ion = 'y';
  }
  var pieces = ion_type.split('(')
  if (pieces.length > 1) {
    var charge = pieces[1][0];
  } else {
    var charge = 1
  }
  var theory_peaks = calculate_peaks(terminii_aa_masses, charge, ion);
  if (terminii_aa_masses.length-2 != sequence.length) {
    console.log('discrepancy in aa_masses', sequence, terminii_aa_masses, modified_aa_masses);
  }
  var matched = [];
  for (var j=0; j<theory_peaks.length; j++) {
    for (var i=0; i<peaks.length; i++) {
      var intensity = peaks[i][1];
      var mz = peaks[i][0];
      var theory_mz = theory_peaks[j][0];
      var label = theory_peaks[j][1];
      var diff = theory_mz - mz;
      var error = Math.round(diff/theory_mz*1000000);
      if (Math.abs(diff) <= mz_delta) {
        matched.push([theory_mz, intensity, label, diff, error]);
        break;
      }
    }
  }
  return matched;
}

seq = 'MSAFLLTKR';
spectrum = [[622.494, 2192.0],
 [615.666, 1811.0],
 [619.456, 1810.0],
 [628.245, 1772.0],
 [616.511, 1016.0],
 [614.583, 780.0],
 [629.003, 673.0],
 [490.215, 563.0],
 [606.378, 533.0],
 [563.797, 430.0],
 [610.871, 421.0],
 [472.861, 376.0],
 [607.103, 346.0],
 [584.56, 341.0],
 [605.67, 334.0],
 [620.333, 310.0],
 [375.222, 299.0],
 [625.311, 263.0],
 [488.298, 254.0],
 [786.349, 223.0],
 [899.377, 221.0],
 [506.51, 215.0],
 [881.571, 215.0],
 [601.8, 207.0],
 [514.999, 202.0],
 [598.419, 201.0],
 [608.138, 199.0],
 [593.698, 196.0],
 [759.359, 196.0],
 [530.555, 186.0],
 [623.222, 186.0],
 [602.825, 184.0],
 [542.091, 180.0],
 [261.063, 169.0],
 [627.002, 160.0],
 [913.689, 160.0],
 [592.707, 157.0],
 [825.581, 157.0],
 [471.826, 156.0],
 [505.835, 155.0],
 [609.168, 154.0],
 [768.284, 151.0],
 [563.027, 146.0],
 [570.712, 146.0],
 [580.803, 145.0],
 [979.473, 137.0],
 [596.779, 130.0],
 [541.113, 129.0],
 [777.839, 186.0],
 [550.083, 159.0]];

// delta_mass = 0.8
// modified_aa_masses = [];
// aa_masses = parse_terminal_aa_masses(seq, aa_monoisotopic_mass, modified_aa_masses);
// console.log(aa_masses)
// modified_aa_masses = [{ 'i':3, 'mass':100.0 }];
// aa_masses = parse_terminal_aa_masses(seq, aa_monoisotopic_mass, modified_aa_masses);
// console.log('modified')
// console.log(modified_aa_masses)
// console.log(aa_masses);
// bion_peaks = calculate_peaks(aa_masses, 1, 'b');
// yion_peaks = calculate_peaks(aa_masses, 1, 'y');
// console.log(bion_peaks.concat(yion_peaks));
// matched = map_matched_ions('y', seq, spectrum, delta_mass, modified_aa_masses, aa_monoisotopic_mass)
// console.log(matched);

