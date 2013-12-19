


function DataController(data) {
  // this is a bit circular but useful
  this.data = data;
  
  this.init = function() {
    if (this.data.mask_labels.length > 0) {
      this.data.mask = parseFloat(this.data.mask_labels[0]);
    } else {
      this.data.mask = 1.0;
    }
    console.log('this.data.mask ' + this.data.mask);
    this.data.canvas_font = "10px 'Andale Mono'";
    this.data.select_bg_color = '#CFC';
    this.data.bg_color = '#F9F9F9';
    this.data.text_color = "#999";
    this.data.selected_seqid = null;
    this.data.delta_mz = 0.5
    this.data.start = true;
    this.data.mass_units = ['m/z(Da)', 'error(ppm)'];
    this.data.mass_unit = this.data.mass_units[0];
    this.data.ion_types = {
       "b(3+)": false,
       "b(2+)": false,
       "b": true,
       "y": true,
       "y(2+)": false,
       "y(3+)": false,
    };
  }

  this.count_peptides = function() {
    for (var seqid in this.data.proteins) {
      var protein = this.data.proteins[seqid];
      protein.attr.n_peptide = 0;
      protein.attr.n_unique_peptide = 0;
      protein.attr.n_slice_populated = 0;
      var sources = protein.sources;
      for (var j=0; j<sources.length; j++) {
        var peptides = sources[j].peptides;
        var n_peptide_in_slice = 0;
        for (var i=0; i<peptides.length; i++) {
          var peptide = peptides[i];
          if (this.data.mask >= peptide.mask) {
            if (peptides[i].attr.is_unique) {
              protein.attr.n_unique_peptide += 1;
            }
            protein.attr.n_peptide += 1;
            n_peptide_in_slice += 1;
          }
        }
        if (n_peptide_in_slice > 0) {
          protein.attr.n_slice_populated += 1;
        }
      }
    }    
  }

  this.check_data = function() {
    for (var seqid in this.data.proteins) {
      if (this.data.selected_seqid == null) {
        this.data.selected_seqid = seqid;
      }
      var protein = this.data.proteins[seqid];
      protein.length = protein.sequence.length;
      protein.i_res_view = 0;
      protein.i_source_view = 0;
      protein.i_peptide_selected = 0;
      protein.i_source_selected = 0;
      for (var i=0; i<protein.sources.length; i++) {
        var peptides = protein.sources[i].peptides;
        for (var j=0; j<peptides.length; j++) {
          var peptide = peptides[j];
          if (('modifications' in peptide.attr) && (peptide.attr.modifications.length == 0)) {
            delete peptide.attr.modifications;
          }
          if (!('j' in peptide)) {
            peptide.j = peptide.i + peptide.sequence.length;
          }
        }
      }
      for (var i=0; i<protein.sources.length; i++) {
        var peptides = protein.sources[i].peptides;
        if (peptides.length > 0) {
          protein.i_source_view = i;
          protein.i_source_selected = i;
          protein.i_res_view = peptides[0]['i'];
          break;
        }     
      }
    }
  }

  this.get_current_protein = function() {
    var seqid = this.data.selected_seqid;
    if (seqid == null) {
      for (key in this.data.proteins) {
        seqid = key;
        break;
      }
    }
    return this.data.proteins[seqid];
  }

  this.get_selected_peptide = function() {
    var protein = this.get_current_protein();
    var i_source = protein.i_source_selected;
    var peptides = protein.sources[i_source].peptides;
    if (peptides.length == 0) {
      return null;
    }
    var i_peptide = protein.i_peptide_selected;
    return peptides[i_peptide];
  }

  this.pick_source_view = function(protein, i_source) {
    protein.i_source_view = i_source;
  }

  this.pick_peptide = function(protein, i_source, i_peptide) {
    protein.i_source_selected = i_source;
    protein.i_peptide_selected = i_peptide;
  }

  this.pick_peptide_callback = function(protein, i_source, i_peptide) {
    var _this = this;
    return function () {
      _this.pick_peptide(protein, i_source, i_peptide);
      _this.data.observer();
      return false;
    }
  }

  this.pick_protein = function(seqid) {
    this.data.selected_seqid = seqid;
  }

  this.pick_protein_callback = function(seqid) {
    var _this = this;
    return function() {
      _this.pick_protein(seqid);
      _this.data.observer();
    }
  }

  this.pick_slice = function(i_source) {
    var protein = this.get_current_protein();
    var peptides = protein.sources[i_source].peptides;
    protein.i_source_view = i_source;
    if (peptides.length > 0) {
      var i_source_old = protein.i_source_selected;
      var i_peptide_old = protein.i_peptide_selected;
      var i_res_old = protein.sources[i_source_old].peptides[i_peptide_old].i;
      var best_i_peptide = 0;
      var best_i_res_diff = Math.abs(peptides[best_i_peptide].i - i_res_old)
      for (var i_peptide=0; i_peptide<peptides.length; i_peptide++) {
        var i_res_diff = Math.abs(peptides[i_peptide].i - i_res_old);
        if (i_res_diff < best_i_res_diff) {
          best_i_peptide = i_peptide;
          best_i_res_diff = i_res_diff;
        }
      }
      protein.i_peptide_selected = best_i_peptide;
      protein.i_source_selected = i_source;
    }
  }

  this.toggle_zoom = function() {
    this.data.zoom = !this.data.zoom;
  }

  this.pick_peak = function(i_peak) {
    this.data.i_peak = i_peak;
  }

  this.pick_peak_callback = function(i_peak) {
    var _this = this;
    return function() {
      if (i_peak != _this.data.i_peak) {
        _this.data.zoom = true;
      } else {
        _this.data.controller.toggle_zoom();
      }
      _this.pick_peak(i_peak);
      _this.data.observer();
      return false;
    }
  }

  this.toggle_ion = function(ion_type) {
    this.data.ion_types[ion_type] = !this.data.ion_types[ion_type];
    this.data.zoom = false;
    var peptide = this.get_selected_peptide();
    delete peptide.labeled_peaks;
  }

  this.toggle_ion_callback = function(ion_type) {
    var _this = this;
    return function() {
      _this.toggle_ion(ion_type);
      _this.data.observer();
      return false;
    }
  }

  this.get_labeled_spectrum = function() {
    var peptide = this.get_selected_peptide();
    if (!('labeled_peaks' in peptide)) {
      console.log('label time');
      peptide.labeled_peaks = peptide.spectrum.slice(0);
      for (ion_type in this.data.ion_types) {
        if (this.data.ion_types[ion_type]) {
          var modifications = [];
          if ('modifications' in peptide.attr) {
             modifications = peptide.attr.modifications;
          }
          var matched = map_matched_ions(
            ion_type,
            peptide.sequence, 
            peptide.labeled_peaks, 
            this.data.delta_mz, 
            modifications, 
            aa_monoisotopic_mass);
          peptide.labeled_peaks = peptide.labeled_peaks.concat(matched);
        }
      }
    }
    return peptide.labeled_peaks;
  }

  this.calc_sorted_seqids = function(sort_key, direction) {
    var sort_pairs = [];
    for (var seqid in this.data.proteins) {
      var protein = this.data.proteins[seqid];
      var sort_val = protein.attr[sort_key];
      sort_pairs.push([sort_val, seqid]);
    }
    if (direction == "smallest") {
      if (typeof sort_pairs[0][0] == 'number') {
        sort_pairs.sort(function(a,b) { return a[0]-b[0] });
      } else {
        sort_pairs.sort()
      }
    } else {
      if (typeof sort_pairs[0][0] == 'number') {
        sort_pairs.sort(function(a,b) { return b[0]-a[0] });
      } else {
        sort_pairs.sort();
        sort_pairs.reverse();
      }
    }
    this.data.sorted_seqids = [];
    for (var i=0; i<sort_pairs.length; i++) {
      this.data.sorted_seqids.push(sort_pairs[i][1]);
    }
  }

  this.get_i_protein = function() {
    var seqid = this.data.selected_seqid;
    for (var j=0; j<this.data.sorted_seqids.length; j++) {
      if (seqid == this.data.sorted_seqids[j]) {
        return j;
      }
    }
    return -1;
  }

  this.onkeydown = function(c) {
    if (c == 'N') {
      var i = this.get_i_protein() + 1;
      if (i < this.data.sorted_seqids.length) {
        var seqid = this.data.sorted_seqids[i];
        this.pick_protein_callback(seqid)();
        return false;
      }
    } else if (c == 'P') {
      var i = this.get_i_protein() - 1;
      if (i >= 0) {
        var seqid = this.data.sorted_seqids[i];
        this.pick_protein_callback(seqid)();
        return false;
      }
    } else if (c == 'J') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      var i_peptide = protein.i_peptide_selected + 1;
      var peptides = protein.sources[i_source].peptides;
      if (i_peptide < peptides.length) {
        this.pick_peptide_callback(protein, i_source, i_peptide)();
      }
    } else if (c == 'K') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      var i_peptide = protein.i_peptide_selected - 1;
      if (i_peptide >= 0) {
        this.pick_peptide_callback(protein, i_source, i_peptide)();
      }
    } else if (c == 'D') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_view;
      if (i_source < protein.sources.length-1) {
        this.pick_slice(i_source+1);
        this.data.observer();
        return false;
      }
    } else if (c == 'U') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_view;
      if (i_source > 0) {
        this.pick_slice(i_source-1);
        this.data.observer();
        return false;
      }
    };
    return;
  }

  this.init();
  this.count_peptides();
  this.check_data();
}

