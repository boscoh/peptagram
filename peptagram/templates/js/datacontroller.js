


function DataController(data) {
  // this is a bit circular but useful
  this.data = data;
  
  this.init = function() {
    if (this.data.mask_labels.length > 0) {
      this.data.mask = parseFloat(this.data.mask_labels[0]);
    } else {
      this.data.mask = 0.0;
    }
    var _data = this.data
    this.data.match_filter = function(match) {
      return match.mask >= _data.mask;
    }
    this.data.canvas_font = "10px 'Andale Mono'";
    this.data.select_bg_color = '#CFC';
    this.data.bg_color = '#F9F9F9';
    this.data.text_color = "#999";
    this.data.selected_seqid = null;
    this.data.delta_mz = 0.5
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

  this.check_data = function() {
    for (var seqid in this.data.proteins) {
      if (this.data.selected_seqid == null) {
        this.data.selected_seqid = seqid;
      }
      var protein = this.data.proteins[seqid];
      protein.length = protein.sequence.length;
      protein.i_res_view = 0;
      protein.i_match_selected = 0;
      protein.i_source_selected = 0;
      for (var i=0; i<protein.sources.length; i++) {
        var matches = protein.sources[i].matches;
        for (var j=0; j<matches.length; j++) {
          var match = matches[j];
          if (('modifications' in match) && (match.modifications.length == 0)) {
            delete match.modifications;
          }
          if (!('j' in match)) {
            match.j = match.i + match.sequence.length;
          }
        }
      }
      for (var i=0; i<protein.sources.length; i++) {
        var matches = protein.sources[i].matches;
        if (matches.length > 0) {
          protein.i_source_selected = i;
          protein.i_res_view = matches[0]['i'];
          break;
        }     
      }
    }
  }

  this.check_location_hash = function() {
    this.data.start = true;
    var hash = window.location.hash.substr(1);
    var params = hash.split('&');

    var pieces = params[0].split('=');
    if (pieces.length < 2) {
      return;
    }
    var seqid = pieces[1];
    if (!(seqid in this.data.proteins)) {
      return;
    }

    this.data.start = false;
    this.data.selected_seqid = seqid;
    var protein = this.data.proteins[seqid];

    var i_source = params[1].split('=')[1];
    var sources = protein.sources;
    if (i_source >= sources.length) {
      i_source = sources.length-1;
    }
    protein.i_source_selected = i_source;

    var i_match = params[2].split('=')[1];
    matches = sources[i_source];
    if (i_match >= matches.length) {
      i_match = matches.length - 1;
    }
    protein.i_match_selected = i_match;
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

  this.get_selected_match = function() {
    var protein = this.get_current_protein();
    var i_source = protein.i_source_selected;
    var matches = protein.sources[i_source].matches;
    var i_match = protein.i_match_selected;
    if ((i_match < 0) || (matches.length == 0)) {
      return null;
    }
    return matches[i_match];
  }

  this.set_location_hash = function() {
    var hash = '#';
    var protein = this.get_current_protein()
    hash += 'protein=' + this.data.selected_seqid;
    hash += '&source=' + protein.i_source_selected;
    hash += '&match=' + protein.i_match_selected;
    window.location.hash = hash;
  }

  this.pick_source_view = function(protein, i_source) {
    protein.i_source_selected = i_source;
    protein.i_match_selected = -1;
    this.set_location_hash();
  }

  this.pick_match = function(protein, i_source, i_match) {
    protein.i_source_selected = i_source;
    protein.i_match_selected = i_match;
    this.set_location_hash();
  }

  this.pick_match_callback = function(protein, i_source, i_match) {
    var _this = this;
    return function () {
      _this.pick_match(protein, i_source, i_match);
      _this.data.observer();
      return false;
    }
  }

  this.pick_protein = function(seqid) {
    this.data.selected_seqid = seqid;
    this.set_location_hash();
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
    var i_source_old = protein.i_source_selected;
    var i_match_old = protein.i_match_selected;
    var new_matches = protein.sources[i_source].matches;
    var old_matches = protein.sources[i_source_old].matches;
    if ((old_matches.length == 0) ||
        (new_matches.length == 0) ||
        (protein.i_match_selected == -1)) {
      protein.i_match_selected = -1;
    } else {
      var i_res_old = old_matches[i_match_old].i;
      var best_i_match = 0;
      var best_i_res_diff = Math.abs(new_matches[best_i_match].i - i_res_old)
      for (var i_match=0; i_match<new_matches.length; i_match++) {
        var i_res_diff = Math.abs(new_matches[i_match].i - i_res_old);
        if (i_res_diff < best_i_res_diff) {
          best_i_match = i_match;
          best_i_res_diff = i_res_diff;
        }
      }
      protein.i_match_selected = best_i_match;
    }
    protein.i_source_selected = i_source;
    this.set_location_hash();
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
    var match = this.get_selected_match();
    delete match.labeled_peaks;
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
    var match = this.get_selected_match();
    if (!('labeled_peaks' in match)) {
      match.labeled_peaks = match.spectrum.slice(0);
      for (ion_type in this.data.ion_types) {
        if (this.data.ion_types[ion_type]) {
          var modifications = [];
          if ('modifications' in match) {
             modifications = match.modifications;
          }
          var matched = map_matched_ions(
            ion_type,
            match.sequence, 
            match.labeled_peaks, 
            this.data.delta_mz, 
            modifications, 
            aa_monoisotopic_mass);
          match.labeled_peaks = match.labeled_peaks.concat(matched);
        }
      }
    }
    return match.labeled_peaks;
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
    if (c == 'A') {
      var i = this.get_i_protein() + 1;
      if (i < this.data.sorted_seqids.length) {
        var seqid = this.data.sorted_seqids[i];
        this.pick_protein_callback(seqid)();
        return false;
      }
    } else if (c == 'W') {
      var i = this.get_i_protein() - 1;
      if (i >= 0) {
        var seqid = this.data.sorted_seqids[i];
        this.pick_protein_callback(seqid)();
        return false;
      }
    } else if (c == 'D') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      var i_match = protein.i_match_selected + 1;
      var matches = protein.sources[i_source].matches;
      if (i_match < matches.length) {
        this.pick_match_callback(protein, i_source, i_match)();
      }
    } else if (c == 'R') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      var i_match = protein.i_match_selected - 1;
      if (i_match >= 0) {
        this.pick_match_callback(protein, i_source, i_match)();
      }
    } else if (c == 'S') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      if (i_source < protein.sources.length-1) {
        this.pick_slice(i_source+1);
        this.data.observer();
        return false;
      }
    } else if (c == 'E') {
      var protein = this.get_current_protein();
      var i_source = protein.i_source_selected;
      if (i_source > 0) {
        this.pick_slice(i_source-1);
        this.data.observer();
        return false;
      }
    };
    return;
  }

  this.init();
  this.check_data();
  this.check_location_hash();
}

