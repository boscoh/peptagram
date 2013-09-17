
//
// Pepto.js a javascript app for displaying multiple 
// Label-Free protein identification using MS/MS analysis
// (c) 2013 Bosco Ho


function attr_empty(dict, key) {
  if (!(key in dict)) {
    return true;
  }
  if (dict[key].length == 0) {
    return true;
  }
  return false;
}


function get_current_protein(data) {
  var seqid = data.selected_seqid;
  if (seqid == null) {
    for (key in data.proteins) {
      seqid = key;
      break;
    }
  }
  return data.proteins[seqid];
}


function get_selected_peptide(data) {
  var protein = get_current_protein(this.data);
  var i_source = protein.i_source_selected;
  var peptides = protein.sources[i_source].peptides;
  if (peptides.length == 0) {
    return null;
  }
  var i_peptide = protein.i_peptide_selected;
  return peptides[i_peptide];
}


function Controller(data) {
  this.data = data;
  
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
    var protein = get_current_protein(this.data);
    var n_source = protein.sources.length;
    var i_source_old = protein.i_source_selected;
    var i_peptide_old = protein.i_peptide_selected;
    var i_res_old = protein.sources[i_source_old].peptides[i_peptide_old].i;
    protein.i_source_view = i_source;
    protein.i_source_selected = i_source;
    var peptides = protein.sources[i_source].peptides;
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
    var peptide = get_selected_peptide(data);
    delete peptide['peaks'];
  }

  this.toggle_ion_callback = function(ion_type) {
    var _this = this;
    return function() {
      _this.toggle_ion(ion_type);
      _this.data.observer();
      return false;
    }
  }

  this.match_ions = function() {
    var peptide = get_selected_peptide(this.data);
    if (!('peaks' in peptide)) {
      peptide.peaks = peptide.spectrum.slice(0);
      for (ion_type in this.data.ion_types) {
        if (this.data.ion_types[ion_type]) {
          var matched = map_matched_ions(
            ion_type,
            peptide.sequence, peptide.peaks, this.data.delta_mz, 
            peptide.attr.modifications, aa_monoisotopic_mass);
          peptide.peaks = peptide.peaks.concat(matched);
        }
      }
    }
  }
}


function make_td(s) {
  return $("<td>").append(s);
}


function draw_highlight_box(canvas, x, y, w, h, color) {
  canvas.box(x-1, y-1, w+2, h+2, "white", 1);
  canvas.box(x-3, y-3, w+6, h+6, color, 3);
}


function sequence_number_line(n_res) {
  var s = "01        ";
  for (var i=10; i<n_res; i+=10) {
    var num = "" + i;
    var n_space = 10 - num.length;
    for (var j=0; j<n_space; j++) {
      num += " ";
    }
    s = s + num;
  }
  s = s.substring(1, n_res+1);
  return s;
}


function sorted_keys(dict) {
  var keys = [];
  for (key in dict) {
    keys.push(key);
  }
  keys.sort();
  return keys;
}


function dict_html(dict) {
  var s = '';
  var keys = sorted_keys(dict);
  for (var i=0; i<keys.length; i++) {
    var key = keys[i];
    var val = dict[key];
    if (key == 'modifications') {
      val = '';
      var modifications = dict[key];
      for (var j=0; j<modifications.length; j++) {
        if (j > 0) {
          val += ', '
        }
        var res_num = modifications[j].i + 1;
        val += res_num + "-" + modifications[j].mass;
      }  
    }
    s += key + ": " + val + "<br>";
  }
  return s;
}


function set_outer_height(div, height) {
  var margin = div.outerHeight(true) - div.innerHeight();
  margin += parseInt(div.css('padding-top'));
  margin += parseInt(div.css('padding-bottom'));
  div.height(height - margin);
}


function set_outer_width(div, width) {
  var margin = div.outerWidth(true) - div.innerWidth();
  margin += parseInt(div.css('padding-left'));
  margin += parseInt(div.css('padding-right'));
  div.width(width - margin);
}


function get_content_width(div) {
  var width = div.innerWidth();
  width -= parseInt(div.css('padding-left'));
  width -= parseInt(div.css('padding-right'));
  return width;
}


function get_content_height(div) {
  var height = div.innerHeight();
  height -= parseInt(div.css('padding-top'));
  height -= parseInt(div.css('padding-bottom'));
  return height;
}


function get_bottom(div) {
  return div.position().top + div.outerHeight(true);
}


// Converts intensities [-1, 1] into an RGB color string
function ColorIntensityPalette() {
  this.neutral = [200, 200, 200];
  this.positive = [242, 90, 90];
  this.negative = [90, 144, 120];
  this.positive_extremum = [255, 0, 0];
  this.negative_extremum = [0, 144, 93];

  this.rgb = function(color) {
      rgb = "rgb(" + 
            color[0] + "," +
            color[1] + "," +
            color[2] + ")";
      return rgb;
  }

  this.str = function(intensity) {
    if (intensity > 1) {
      return this.rgb(this.positive_extremum);
    } else if (intensity < -1) {
      return this.rgb(this.negative_extremum);
    }

    if (intensity > 0) {
      var top_color = this.positive;
    } else {
      var top_color = this.negative;
      intensity = -intensity;
    }

    s = "rgb(";
    for (var i=0; i<3; i++) {
      diff = top_color[i] - this.neutral[i];
      s += Math.floor(this.neutral[i] + diff*intensity);
      if (i < 2) {
        s += ","
      }
    }
    s += ")"
    return s;
  }
}


function ColorBarWidget(canvas, names, font) {
  this.canvas = canvas;
  this.font = font;
  this.names = names;
  this.palette = new ColorIntensityPalette();

  this.x = 0;
  this.y = 0;

  this.name_width = 0;
  this.canvas.draw_context.font = this.font;
  for (var i=0; i<this.names.length; i++) {
    var name = this.names[i];
    var test_name_width = this.canvas.draw_context.measureText(name).width;
    if (test_name_width > this.name_width) {
      this.name_width = test_name_width;
    }
  }
  this.bar_width = 10;
  this.width = this.bar_width + 2 + this.name_width;
  this.height = 80;

  this.draw = function() {
    var height = this.height - 10;
    var half_height = height/2;
    var y_min = this.y + 5;
    var x2 = this.x + this.bar_width;
    var y_max = height;
    for (var i=0; i<height; i++) {
      f = (half_height - i)/half_height;
      y = y_min + i;
      this.canvas.line(this.x, y, x2, y, this.palette.str(f), 1);
    }
    var canvas = this.canvas;
    var y = y_min - 5;
    this.canvas.line(this.x, y, x2, y, this.palette.str(2), 10);
    var y = y_min + height + 5;
    this.canvas.line(this.x, y, x2, y, this.palette.str(-2), 10);
    var write_name = function(txt, y) {
      canvas.text(txt, x2 + 2, y, this.font, "#999", "left");
    }
    write_name(this.names[0], this.y);
    write_name(this.names[1], this.y + 0.5*this.height);
    write_name(this.names[2], this.y + this.height);
  }
}


// ProteinBarWidget draws a simple distribution of peptides in a
// protein using data.mask to mask certain peptides
function ProteinBarWidget(canvas, data, seqid) {
  this.canvas = canvas;
  this.seqid = seqid;
  this.data = data;
  this.protein = this.data.proteins[seqid];
  this.color = "#999";
  this.x = 0;
  this.y = 0;
  this.width = this.canvas.canvas_dom.width;
  this.height = this.canvas.canvas_dom.height;

  this.x_from_i = function(i) {
    return i/this.protein.length*this.width + this.x;
  }

  this.draw = function() {
    if (this.seqid == this.data.selected_seqid) {
      var bg_color = this.data.select_bg_color;
    } else {
      var bg_color = '#EEE';
    }
    this.canvas.solid_box(this.x, this.y, this.width, this.height, bg_color);
    this.canvas.solid_box(this.x, this.y + this.height/2, this.width, 1, '999');
    for (var j=0; j<this.protein.sources.length; j++) {
      var peptides = this.protein.sources[j].peptides;
      for (var i=0; i<peptides.length; i++) {
        var peptide = peptides[i];
        if (this.data.mask <= peptide.mask) {
          var x = this.x_from_i(peptide.i);
          var w = this.x_from_i(peptide.j) - x;
          this.canvas.solid_box(x, this.y, w, this.height, this.color);
        }
      }
    }
  }
}


function ProteinList(control_div, column1_div, data) {
  this.main_div = column1_div
  this.control_div = control_div;
  this.data = data;
  this.protein_widgets = [];
  this.protein_divs = [];
  this.key;
  this.sorting_msg_div;
  this.attr_select;

  this.build_controls = function() {
    // selectors for sorting protein attributes
    this.control_div.append("Sort by:");
    var keys = sorted_keys(get_current_protein(this.data).attr);
    this.attr_select = $("<select>");
    for (i=0; i<keys.length; i++) {
      var key = keys[i];
      var option = $('<option>');
      option.attr('value', key);
      option.text(key);
      if (key == 'n_peptide') {
        option.attr('selected', true);
      }
      this.attr_select.append(option);
    }
    this.control_div.append(this.attr_select);
    this.control_div.append("<br>");

    // radio buttons for Largest/Smallest
    var directions = ['largest', 'smallest'];
    for (var i=0, direction; direction=directions[i]; i++) {
      var button = $("<input>");
      button.attr('name', 'direction');
      button.attr('type', 'radio');
      button.attr('value', direction);
      if (i == 0) {
        button.attr('checked', true);
      }
      this.control_div.append(button)
      this.control_div.append(direction + '&nbsp;');
    }
    this.control_div.append('<br>');

    // radio buttons for masking
    if (this.data['mask_labels'].length > 0) {
      this.control_div.append('FPE: ');
      for (var i=0; i<this.data['mask_labels'].length; i++) {
        var button = $("<input>");
        button.attr('name', 'mask');
        button.attr('type', 'radio');
        button.attr('value', '' + i);
        if (i == 0) {
          button.attr('checked', true);
        }
        this.control_div.append(button);
        this.control_div.append(this.data['mask_labels'][i] + ' ');
      }
      this.control_div.append('<br>');
    }
    this.control_div.append('<br>');
  }

  this.calc_sorted_seqids = function() {
    var unsorted_seqids = [];
    this.key = this.attr_select.val();
    for (var seqid in this.data.proteins) {
      var protein = this.data.proteins[seqid];
      var sort_key = protein.attr[this.key];
      unsorted_seqids.push([sort_key, seqid]);
    }
    var direction = $("input[name=direction]:checked").attr('value');
    if (direction == "smallest") {
      if (typeof unsorted_seqids[0][0] == 'number') {
        unsorted_seqids.sort(function(a,b) { return a[0]-b[0] });
      } else {
        unsorted_seqids.sort()
      }
    } else {
      if (typeof unsorted_seqids[0][0] == 'number') {
        unsorted_seqids.sort(function(a,b) { return b[0]-a[0] });
      } else {
        unsorted_seqids.sort();
        unsorted_seqids.reverse();
      }
    }
    this.data.sorted_seqids = [];
    for (var i=0; i<unsorted_seqids.length; i++) {
      this.data.sorted_seqids.push(unsorted_seqids[i][1]);
    }
  }

  this.build_sorting_msg = function() {
    this.sorting_msg_div = $('<div>')
      .text("SORTING: PROTEIN LIST...")
      .addClass('title')
      .addClass('sorting_msg');
    this.control_div.append(this.sorting_msg_div)
    var width = get_content_width(this.control_div);
    set_outer_width(this.sorting_msg_div, width);
  }

  this.get_i_protein = function(seqid) {
    for (var j=0; j<this.data.sorted_seqids.length; j++) {
      if (seqid == this.data.sorted_seqids[j]) {
        return j;
      }
    }
    return -1;
  }

  this.highlight_selected_protein = function() {
    this.i_protein = this.get_i_protein(this.data.selected_seqid);
    if (this.i_old_protein == this.i_protein) {
      return;
    }
    this.data.selected_seqid = this.data.sorted_seqids[this.i_protein];
    if (this.i_old_protein != null) {
      this.protein_divs[this.i_old_protein].css('color', '#999');
      this.protein_widgets[this.i_old_protein].draw();
    }
    this.protein_divs[this.i_protein].css({'color':'#600'});
    this.protein_widgets[this.i_protein].draw();
    this.i_old_protein = this.i_protein;
  }


  this.build_list = function() {
    this.main_div.empty();
    this.main_div.append('<br>');
    this.build_sorting_msg();
    this.calc_sorted_seqids();
    while (this.protein_widgets.length > 0) {
      this.protein_widgets.pop();
    }
    while (this.protein_divs.length > 0) {
      this.protein_divs.pop();
    }

    if (this.data.start) {
      this.data.selected_seqid = this.data.sorted_seqids[0];
      this.data.start = false;
      this.data.controller.pick_protein(this.data.selected_seqid);
    }
    
    // Async call-back function to build protein-list
    // incrementally so the page doesn't freeze
    var ms_wait = 20;
    var i_seqid = 0;
    var n_seqid_batch = 20;
    var _this = this;
    this.i_old_protein = null;
    var progressive_build_list = function() {
      var j_seqid = Math.min(_this.data.sorted_seqids.length, i_seqid + n_seqid_batch);
      for (; i_seqid<j_seqid; i_seqid++) {
        var seqid = _this.data.sorted_seqids[i_seqid];
        var protein = _this.data.proteins[seqid];
        var val = protein.attr[_this.key];
        var description = "" + (i_seqid+1) + ":";
        description += "[" + val + "] ";
        description += seqid + ": ";
        description += protein.description;

        var column2_div = $("<div>")
        column2_div.append($("<div>").text(description));
        column2_div.click(_this.data.controller.pick_protein_callback(seqid));
        _this.protein_divs.push(column2_div);
        _this.main_div.append(column2_div);

        var canvas_div = $("<div>");
        canvas_div.css('width', _this.main_div.width() - 10);
        canvas_div.css('height', 12);
        canvas_div.css('margin', '5px 0');
        column2_div.append(canvas_div);

        var canvas = new CanvasWidget(canvas_div, "#EFEFEF");
        var protein_bar_widget = new ProteinBarWidget(canvas, data, seqid);
        _this.protein_widgets.push(protein_bar_widget);

        canvas.push(protein_bar_widget);
        canvas.draw();
      }
      
      if (j_seqid < _this.data.sorted_seqids.length) {
        setTimeout(progressive_build_list, ms_wait);
      } else {
        _this.main_div.append('<br><br>')
        _this.sorting_msg_div.remove();
        _this.highlight_selected_protein();
      }
    }
    progressive_build_list();
  }

  this.redraw_mask = function() {
    var val = $("input[name=mask]:checked").attr('value');
    this.data.mask = parseInt(val);
    this.build_list();
    this.data.observer();
  }

  this.register_callbacks = function() {
    var _this = this;
    this.attr_select.change(
        function() { _this.build_list(); });
    $("input[name=direction]").change(
        function() { _this.build_list(); });
    $("input[name=mask]").change(
        function() { _this.redraw_mask(); });
  }

  this.build_controls();
  this.build_list();
  this.register_callbacks();
}


function build_protein_info_panel(data, div) {
  div.empty();
  var protein = get_current_protein(data);
  div.append(protein['description']);
  div.append("<br>----<br>");
  div.append(dict_html(protein['attr']));
}


function PeptographWidget(canvas, data, color_bar) {
  this.canvas = canvas;
  this.seqid = null;
  this.data = data;
  this.color_bar = color_bar;
  this.pressed = false;
  this.down = false;
  this.up = false; 
  this.x = 0;
  this.y = 0;
  this.protein = get_current_protein(this.data);

  this.x_from_i = function(i) {
    return i/this.protein.length*this.draw_width + this.x;
  }

  this.y_from_j = function(j) {
    return j/this.n_source*this.draw_height + this.y;
  }

  this.i_from_x = function(x) {
    return Math.floor((x - this.x)/this.draw_width*this.protein.length);
  }

  this.j_from_y = function(y) {
    var j = Math.floor((y - this.y)/this.draw_height*this.n_source);
    if (j >= this.n_source) {
      j = this.n_source - 1;
    }
    if (j < 0) {
      j = 0;
    }
    return j;
  }

  this.get_diff_width = function(i0, i1) {
    return Math.abs(this.x_from_i(i1) - this.x_from_i(i0));
  }
  
  this.get_diff_height = function(i0, i1) {
    return Math.abs(this.y_from_j(i1) - this.y_from_j(i0));
  }
  
  this.draw = function() {
    this.protein = get_current_protein(this.data);

    var sources = this.protein.sources;
    this.n_source = sources.length;

    var name_width = 0;
    this.canvas.draw_context.font = this.data.canvas_font;
    for (var i=0; i<this.data.source_labels.length; i++) {
      var label = this.data.source_labels[i];
      var label_width = this.canvas.draw_context.measureText(label).width;
      if (label_width > name_width) {
        name_width = label_width;
      }
    }
    if (name_width > 0) {
      name_width += 10;
    }

    var padding = 4;
    var freq_width = 50;
    var slice_height = this.get_diff_height(1, 0);
    var freq_height = slice_height/2;

    // fit around the color_bar which is implicitly to the left
    this.x = this.color_bar.x + this.color_bar.width + 10 + padding;
    this.y = padding;
    this.draw_width = this.canvas.canvas_dom.width - this.x - name_width - freq_width;
    this.draw_height = this.canvas.canvas_dom.height - this.y - padding;

    // draw background area
    this.canvas.solid_box(
        this.x, this.y, this.draw_width, this.draw_height, "white");

    // draw slice freq sidebar
    var max_n_pep = 0;
    for (var j=0; j<sources.length; j++) {
      if (sources[j].peptides.length > max_n_pep) {
        max_n_pep = sources[j].peptides.length;
      }
    }
    for (var i_source=0; i_source<sources.length; i_source++) {
      var n_pep = sources[i_source].peptides.length;
      if (n_pep > 0) {
        this.canvas.solid_box(
          this.x + this.draw_width, 
          this.y_from_j(i_source) + freq_height/2, 
          n_pep/max_n_pep*freq_width, 
          freq_height, 
          '#999');
      }
    }

    // draw max slice freq num
    this.canvas.text(
        '' + max_n_pep, 
        this.x + this.draw_width + freq_width - 2,
        this.y_from_j(0) + slice_height/2, 
        this.data.canvas_font, 
        "#CCC", 
        "right");

    // draw selected row background
    this.canvas.solid_box(
        this.x, 
        this.y_from_j(this.protein.i_source_view),
        this.draw_width, 
        slice_height, 
        this.data.select_bg_color);

    // draw labels
    if ('source_labels' in this.data) {
      for (var j=0; j<sources.length; j++) {
        if (j < this.data.source_labels.length) {
          this.canvas.text(
              this.data.source_labels[j], 
              this.x + this.draw_width + freq_width + name_width,
              this.y_from_j(j) + slice_height/2, 
              this.data.canvas_font, 
              "#999", 
              "right");
        }
      }
    }

    // draw peptides
    for (var j=0; j<sources.length; j++) {
      var peptides = sources[j].peptides;
      for (var i=0; i<peptides.length; i++) {
        var peptide = peptides[i];
        if (this.data.mask <= peptide.mask) {
          if (peptide.intensity === "") {
            var color = "lightyellow";
          } else {
            var color = this.color_bar.palette.str(peptide.intensity);
          }
          this.canvas.solid_box(
              this.x_from_i(peptide.i), 
              this.y_from_j(j), 
              this.get_diff_width(peptide.i, peptide.j), 
              slice_height, 
              color);
        }
      }
    }

    //  highlight sequence_view area
    draw_highlight_box(
        this.canvas, 
        this.x_from_i(this.protein.i_res_view), 
        this.y_from_j(this.protein.i_source_view), 
        this.get_diff_width(0, this.data.n_res_in_view), 
        slice_height,
        "#DA3");

    //  highlight selected peptide
    var i_source_selected = this.protein.i_source_selected;
    var peptides = sources[i_source_selected].peptides;
    if (peptides.length > 0) {
      peptide = peptides[this.protein.i_peptide_selected];
      draw_highlight_box(
          this.canvas, 
          this.x_from_i(peptide.i),
          this.y_from_j(i_source_selected),
          this.get_diff_width(peptide.i, peptide.j),
          slice_height,
          "#3F3");
    }
  }

  this.drag = function(x, y) {
    if (this.pressed === false) {
      return;
    }
    // set view params
    var i_res = this.i_from_x(x);
    var i_res_view = i_res - this.data.n_res_in_view/2;
    if (i_res_view < 0) {
      i_res_view = 0;
    }
    if ((i_res_view + this.data.n_res_in_view) > this.protein.length) {
      i_res_view = this.protein.length - this.data.n_res_in_view;
    }
    this.protein.i_res_view = i_res_view;
    var i_source = this.j_from_y(y);
    this.data.controller.pick_source_view(this.protein, i_source);

    // if click over peptide: select it
    var source = this.protein.sources[i_source];
    for (var i=0; i<source.peptides.length; i++) {
      var peptide = source.peptides[i];
      if ((peptide.i<=i_res) && (i_res<peptide.j)) {
        this.data.controller.pick_peptide(this.protein, i_source, i);
        break;
      }
    }

    this.data.observer();
  }
}


var SequenceView = function(div, data) {
  this.data = data;
  this.seqid = null;
  this.div = div;
  this.dom = this.div[0];

  this.char_width = function() {
    return this.dom.scrollWidth/this.protein.length;
  }

  this.get_scroll_i_res = function() {
    var n_res = this.protein.length;
    var scroll_left = this.div.scrollLeft();
    var scroll_width = this.dom.scrollWidth;
    return Math.floor(n_res*scroll_left/scroll_width);
  }

  this.build = function() {
    this.seqid = this.data.selected_seqid;
    this.protein = get_current_protein(this.data);
    this.i_source = this.protein.i_source_view;
    this.i_peptide_selected = this.protein.i_peptide_selected;
   
    var peptides = this.protein.sources[this.i_source].peptides;
    var sequence = this.protein.sequence;
    var n_res = sequence.length;

    var is_peptides_intersect = function(i1, i2) {
      var pep1 = peptides[i1];
      var pep2 = peptides[i2];
      var not_intersect = ((pep1.j <= pep2.i) || (pep2.j <= pep1.i));
      return !not_intersect;
    }

    this.div.empty();

    var sequence_header = sequence_number_line(n_res);
    var pre = $('<pre>')
    this.div.append(pre);

    pre.append(sequence_header);
    pre.append('<br>');

    i_res = 0;
    for (var i_peptide=0; i_peptide<peptides.length; i_peptide++) {
      peptide = peptides[i_peptide];
      if (this.data.mask > peptide.mask) {
        continue;
      }
      if (i_res < peptide.i) {
        pre.append(sequence.slice(i_res, peptide.i));
        i_res = peptide.i;
      }
      var peptide_link = $("<a>");
      if (this.i_peptide_selected < peptides.length) {
        if (is_peptides_intersect(i_peptide, this.i_peptide_selected)) {
          peptide_link.addClass('highlight_peptide');
        }
      }
      peptide_link.attr('href', '#');
      peptide_link.addClass('link');
      peptide_link.text(sequence.slice(i_res, peptide.j));
      peptide_link.click(this.data.controller.pick_peptide_callback(this.protein, this.i_source, i_peptide));
      pre.append(peptide_link);
      i_res = peptide.j;
    }
    if (i_res < n_res) {
      pre.append(sequence.slice(i_res, n_res));
    }
    this.data.n_res_in_view = Math.ceil(this.div.width()/this.char_width());
  }

  this.update_panel = function() {
    if ((this.seqid !== this.data.selected_seqid) ||
        (this.i_source !== this.protein.i_source_view) ||
        (this.i_peptide_selected != this.protein.i_peptide_selected)) {
      this.build();
    }
    var i_res = this.protein.i_res_view;
    var scroll_i_res = this.get_scroll_i_res();
    // only scroll if more than one char different
    if (Math.abs(scroll_i_res - i_res) > 0) {
      var new_scroll_left = Math.floor(i_res*this.char_width())
      this.div.scrollLeft(new_scroll_left);
    }
  }

  this.register_callbacks = function() {
    var _this = this;
    this.div.scroll(function() { 
       _this.protein.i_res_view = _this.get_scroll_i_res();
       _this.data.observer();
       return false 
    });
  }

  this.register_callbacks();
}


function build_peptides_panel(data, div) {
  div.empty();
  var protein = get_current_protein(data);
  var i_source = protein.i_source_selected;
  var peptides = protein.sources[i_source].peptides;
  var table = $('<table>');
  table.css('text-align', 'left');
  div.append(table);

  // heading
  var tr = $('<tr>')
  table.append(tr);
  tr.append($('<td>').text('#'));
  tr.append($('<td>').text('i'));
  tr.append($('<td>').text('seq'));

  for (i_peptide=0; i_peptide<peptides.length; i_peptide++) {
    var peptide = peptides[i_peptide];

    var tr = $('<tr>')
    table.append(tr);

    var td = $('<td>');
    td.append(i_peptide+1);
    tr.append(td);

    var td = $('<td>');
    td.append(peptide.i);
    tr.append(td);

    var td = $('<td>');
    var peptide_link = $("<a>");
    if ((i_peptide == protein.i_peptide_selected)) {
      peptide_link.addClass('highlight_peptide');
    }
    peptide_link.attr('href', '#');
    peptide_link.addClass('link');
    var seq = peptide['sequence'];
    var modifications = peptide.attr.modifications;
    var modified = [];
    for (var i=0; i<seq.length; i++) {
      modified.push(false);
    }
    for (var i=0; i<modifications.length; i++) {
      modified[modifications[i].i] = true;
    }
    for (var i=0; i<seq.length; i++) {
      if (modified[i]) {
        var s = $('<span>').append(seq[i]).addClass('modification');
      } else {
        var s = seq[i];    
      }
      peptide_link.append(s);
    }
    peptide_link.click(data.controller.pick_peptide_callback(protein, i_source, i_peptide));
    td.append(peptide_link);
    tr.append(td);
  }
}


function build_peptide_info(data, div) {
  div.empty();
  var peptide = get_selected_peptide(data);
  div.append(dict_html(peptide['attr']));
  div.append("<br>");
}


function SpectrumWidget(canvas, data) {
  this.canvas = canvas;
  this.data = data;
  this.x = 0;
  this.y = 0;
  this.offset = 30;
  this.pressed = false;
  this.down = false;
  this.up = false; 
  this.data.zoom = false;
  this.data.i_peak = null;

  this.get_optimum_limits = function() {
    this.min_m = 1E6;
    this.max_m = 0;
    this.min_i = 1E6;
    this.max_i = 0;
    for (var j=0; j<this.spectrum.length; j++) {
      var m = this.spectrum[j][0];
      var i = this.spectrum[j][1];
      if (m < this.min_m) {
        this.min_m = m;
      }
      if (m > this.max_m) {
        this.max_m = m;
      }
      if (i < this.min_i) {
        this.min_i = i;
      }
      if (i > this.max_i) {
        this.max_i = i;
      }
    }
    this.diff_m = this.max_m - this.min_m;
    this.diff_i = this.max_i - this.min_i;
    this.min_m -= this.diff_m*0.1;
    this.max_m += this.diff_m*0.1;
    this.diff_m = this.max_m - this.min_m;
    this.diff_i = this.max_i - this.min_i;
  }

  this.get_zoom_limits = function() {
    var peak = this.spectrum[this.data.i_peak];
    var m = peak[0];
    var delta_m = 20;
    this.min_m = m - delta_m;
    this.max_m = m + delta_m;
    this.min_i = 0;
    this.max_i =  peak[1];
    this.diff_m = this.max_m - this.min_m;
    this.diff_i = this.max_i - this.min_i;
  }

  this.x_from_m = function(m) {
    return this.x + this.offset
        + (m-this.min_m)/this.diff_m*(this.draw_width);
  }

  this.m_from_x = function(x) {
    return (x - this.x - this.offset) / 
           this.draw_width*this.diff_m + 
           this.min_m;
  }

  this.y_from_i = function(i) {
    if (i > this.max_i) {
      i = this.max_i;
    }
    return this.y + this.offset + this.draw_height 
        - (i-this.min_i)/this.diff_i*(this.draw_height);
  }

  this.draw = function() {
    this.width = this.canvas.canvas_dom.width;
    this.height = this.canvas.canvas_dom.height;
    this.draw_width = this.width - 2*this.offset;
    this.draw_height = this.height - 2*this.offset;
    this.data.delta_mz = 0.5

    var peptide = get_selected_peptide(this.data);
    if (attr_empty(peptide, 'spectrum')) {
      this.canvas.div.css('display', 'none');
      return;
    } 
    this.canvas.div.css('display', 'block');

    this.data.controller.match_ions();
    this.spectrum = peptide.peaks;
    if (this.data.i_peak >= this.spectrum.length) {
      this.data.i_peak = this.spectrum.length-1;
    }

    canvas.solid_box(this.x, this.y, this.width, this.height, "#EFE");

    if (this.data.zoom == false) {
      this.get_optimum_limits();
      var title = 'CLICK TO ZOOM';
    } else {
      this.get_zoom_limits();
      var peak = this.spectrum[this.data.i_peak];
      var title = 'CLICK TO ZOOM OUT';
    }

    // draw title
    canvas.text(
      title, 
      this.x + this.width/2,
      this.y + 8,
      this.data.canvas_font, '#999', 'center')

    for (var i=0; i<this.spectrum.length; i++) {
      var x1 = this.x_from_m(this.spectrum[i][0]);
      if ((x1 < this.x + this.offset) ||
          (x1 > this.x + this.offset + this.draw_width)) {
        continue;
      }
      var y1 = this.y_from_i(this.spectrum[i][1]);
      var y2 = this.y + this.offset + this.draw_height;
      if (y2 == y1) {
        y1 = y2 - 1;
      }
      if (this.spectrum[i].length > 2) {
        var label = this.spectrum[i][2];
      } else {
        var label = '';
      }
      if (label) {
        if (label[0] == "y") {
          color = "orange";
        } else if (label[0] == "b") {
          color = "purple";
        } 
      } else {
        color = "rgb(200, 200, 200)";
      }
      if (this.data.zoom) {
        canvas.line(x1, y1, x1, y2, color, 3);
      } else {
        canvas.line(x1, y1, x1, y2, color, 1);
      }
      if (label) {
        canvas.text(label, x1, y1 - 8, this.data.canvas_font, color);
      }
    }

    // draw zoom help line
    if (!this.data.zoom) {
      var i = this.data.i_peak;
      if (i) {
        var x1 = this.x_from_m(this.spectrum[i][0]);
         if ((x1 >= this.x + this.offset) &&
            (x1 <= this.x + this.offset + this.draw_width)) {
          var y1 = this.y + this.offset;
          var y2 = this.y + this.offset + this.draw_height;
          canvas.line(x1, y1, x1, y2, '#AFA', 2);
        }
      }
    }

    // draw frame
    canvas.line(
        this.x + this.offset,
        this.y + this.offset,
        this.x + this.offset,
        this.y + this.offset + this.draw_height,
        "#EEE", 2);
    canvas.line(
        this.x + this.offset,
        this.y + this.offset + this.draw_height,
        this.x + this.offset + this.draw_width,
        this.y + this.offset + this.draw_height,
        "#EEE", 2);

    // draw min, max values of mass/charage
    canvas.text(
      this.min_m.toFixed(0), 
      this.x + this.offset, 
      this.y + this.offset + this.draw_height + 10, 
      this.data.canvas_font, '#999', 'center')
    canvas.text(
      this.max_m.toFixed(0), 
      this.x + this.offset + this.draw_width, 
      this.y + this.offset + this.draw_height + 10, 
      this.data.canvas_font, '#999', 'center')
    if (this.data.zoom) {
      var label = '' + peak[0].toFixed(0);
      canvas.text(
        label, 
        this.x + this.offset + this.draw_width/2, 
        this.y + this.offset + this.draw_height + 10, 
        this.data.canvas_font, '#999', 'center')
    }

    // draw intensity max label
    canvas.text(
      this.max_i.toFixed(0), 
      this.x + 2,
      this.y + this.offset - 8,
      this.data.canvas_font, '#999', 'left')
  }

  this.get_i_peak = function(x) {
    var i_peak_closest = 0;
    var diff_closest = 100;
    var m = this.m_from_x(x);
    for (var j=0; j<this.spectrum.length; j++) {
      var peak = this.spectrum[j];
      var diff = Math.abs(peak[0] - m);
      if (diff < diff_closest) {
        diff_closest = diff;
        i_peak_closest = j;
      }
    }
    return i_peak_closest;
  }


  this.drag = function(x, y) {
    var peptide = get_selected_peptide(this.data);
    if (attr_empty(peptide, 'spectrum')) {
      this.canvas.div.css('display', 'none');
      return;
    }
    if (!this.data.zoom) {
      this.data.controller.pick_peak(this.get_i_peak(x));
    }
    if (this.down) {
      this.data.controller.toggle_zoom();
    }
    this.data.observer();
  }
}


function build_ion_table_div(data, div) {
  div.empty();
  var peptide = get_selected_peptide(data);
  if (attr_empty(peptide, 'spectrum')) {
    return;
  }
  var spectrum = peptide.peaks;

  div.append('<br>')
  for (var i=0, mass_unit; mass_unit=data.mass_units[i]; i++) {
    var button = $("<input>");
    button.attr('name', 'mass_unit');
    button.attr('type', 'radio');
    button.attr('value', mass_unit);
    if (data.mass_unit == mass_unit) {
      button.attr('checked', true);
    };
    div.append(button)
    div.append(mass_unit + '&nbsp;');
  }
  $("input[name=mass_unit]").change(function() { 
      var mass_unit = $("input[name=mass_unit]:checked").attr('value');
      data.mass_unit = mass_unit;
      data.observer(); 
  });

  function make_ion_heading(ion_type) {
    var s = ion_type;
    if (s.length == 1) {
      s = '&nbsp;' + s + '&nbsp;';
    }
    var a = $('<a>').html(s).addClass('link').attr('href','')
    var td = make_td(a);
    a.click(data.controller.toggle_ion_callback(ion_type));
    if (data.ion_types[ion_type]) {
      a.addClass('highlight_peptide');
    }
    return td;
  }

  // build table header
  div.append('<hr>')
  var table = $("<table>");
  var tr = $("<tr>");
  tr.append(make_ion_heading("b(3+)"));
  tr.append(make_ion_heading("b(2+)"));
  tr.append(make_ion_heading("b"));
  tr.append(make_td("C"));
  tr.append(make_ion_heading("y"));
  tr.append(make_ion_heading("y(2+)"));
  tr.append(make_ion_heading("y(3+)"));
  table.append(tr);

  var mass_unit = $("input[name=mass_unit]:checked").attr('value');

  function mass_td(label) {
    var s = '';
    for (var i=0; i<spectrum.length; i++) {
      if (spectrum[i][2] === label) {
        s = $('<a>').addClass('link');
        if (mass_unit == data.mass_units[0]) {
          s.text(spectrum[i][0].toFixed(0));
        } else {
          s.text(spectrum[i][4].toFixed(0));
        }
        s.attr('href', '');
        s.click(this.data.controller.pick_peak_callback(i))
        break;
      }
    }
    return make_td(s);
  }

  // build ion table
  var n = peptide.sequence.length
  for (var i=0; i<n; i++) {
    var tr = $("<tr>");
    var num = i + 1;
    tr.append(mass_td(ion_label(num, 3, "b")));
    tr.append(mass_td(ion_label(num, 2, "b")));
    tr.append(mass_td(ion_label(num, 1, "b")));
    pre = '' + num;
    if (pre.length == 1) {
      pre = '.' + pre;
    }
    post = '' + (peptide.sequence.length - i);
    if (post.length == 1) {
      post = post + '.';
    }
    tr.append(make_td(pre + '.' + peptide.sequence[i] + '.' + post));
    var num = peptide.sequence.length - i;
    tr.append(mass_td(ion_label(num, 1, "y")));
    tr.append(mass_td(ion_label(num, 2, "y")));
    tr.append(mass_td(ion_label(num, 3, "y")));
    table.append(tr);
  }
  div.append(table);
  div.append("<hr>");

  // build extra ions with prosthetic groups
  s = "";
  for (var i=0; i<spectrum.length; i++) {
    if (spectrum[i].length > 2) {
      var word = spectrum[i][2];
    } else {
      word = '';
    }
    if (word.indexOf("-") >= 0) {
      s += spectrum[i][2]; 
      s += "(" + spectrum[i][0].toFixed(0) + "), ";
    }
  }
  if (s) {
    div.append(s);
    div.append("<hr>");
  }
}


function Pepto(data) {
  this.init_data = function(data) {
    this.data = data;
    this.data.mask = 0;
    this.data.canvas_font = this.data.canvas_font;
    this.data.select_bg_color = '#FED';
    this.data.selected_seqid = null;
    this.data.start = true;
    this.data.mass_units = ['mz(Da)', 'error(ppm)'];
    this.data.mass_unit = this.data.mass_units[0];
    this.data.ion_types = {
       "b(3+)": false,
       "b(2+)": false,
       "b": true,
       "y": true,
       "y(2+)": false,
       "y(3+)": false,
    };
    this.data.n_protein = 0;
    this.data.controller = new Controller(this.data);
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
          if (this.data.mask <= peptide.mask) {
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
    for (var seqid in this.data.proteins) {
      if (this.data.selected_seqid == null) {
        this.data.selected_seqid = seqid;
      }
      this.data.n_protein += 1;
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
          if (!('modifications' in peptide.attr)) {
            peptide.attr.modifications = [];
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

  this.init_widgets = function() {
    this.header_div = $('#header');
    this.column1_div = $('#column1');
    this.protein_control_div = $("#protein_control");
    this.protein_list_div = $("#protein_list");
    this.column2_div = $('#column2');
    this.protein_info_div = $('#protein_info')
    this.peptograph_div = $("#peptograph")
    this.sequence_div = $("#sequence");
    this.column3_div = $("#column3");
    this.peptide_list_div = $("#peptide_list");
    this.peptide_info_div = $('#peptide_info')
    this.ion_table_div = $('#ion_table')
    this.spectrum_div = $("#spectrum");
    this.peptide_match_info_div = $("#peptide_match_info");

    this.header_div.text(this.data['title']);

    // Build interactive page objects
    this.protein_list = new ProteinList(
        this.protein_control_div, this.protein_list_div, this.data);

    this.peptograph_canvas = new CanvasWidget(this.peptograph_div, "#EEE");
    var width = this.peptograph_canvas.div.width();
    this.color_bar = new ColorBarWidget(
        this.peptograph_canvas, this.data['color_names'], this.data.canvas_font);
    this.color_bar.x = 0;
    this.color_bar.y = 5;
    var peptograph_widget = new PeptographWidget(
        this.peptograph_canvas, this.data, this.color_bar);
    peptograph_widget.x = 30;
    this.peptograph_canvas.push(this.color_bar);
    this.peptograph_canvas.push(peptograph_widget);

    this.sequence_view = new SequenceView(this.sequence_div, this.data);

    this.spectrum_canvas = new CanvasWidget(this.spectrum_div, "#EFEFEF");
    spectrum_widget = new SpectrumWidget(this.spectrum_canvas, this.data);
    this.spectrum_canvas.push(spectrum_widget);

  }

  this.resize_display = function() {
    var window_width = $(window).width();
    var window_height = $(window).height();

    // set columns at the right height
    var header_height = this.header_div.outerHeight(true);
    this.column1_div.css('top', header_height);
    this.column2_div.css('top', header_height);
    this.column3_div.css('top', header_height);
    var main_height = window_height - header_height;

    // put the protein widget bars in the right place
    var height = this.column1_div.innerHeight();
    var top =  get_bottom(this.protein_control_div);
    set_outer_height(this.protein_list_div, height-top);
    set_outer_height(this.column1_div, main_height);

    // move the central column to its right place
    var protein_list_width = this.column1_div.outerWidth(true);
    this.column2_div.css('left', protein_list_width);

    // set heights and tops of central protein view
    set_outer_height(this.column2_div, main_height);
    var top = get_bottom(this.protein_info_div);
    this.peptograph_div.css('top', top);
    var height = get_content_height(this.column2_div) - 1;
    height += - top - this.sequence_div.outerHeight(true);
    var n_source = get_current_protein(this.data).sources.length;
    if (n_source < 5) {
      var fixed_height = this.color_bar.height*n_source + 10;
      if (fixed_height < height) {
        height = fixed_height;
      }
    }
    this.peptograph_canvas.set_height(height);
    top = get_bottom(this.peptograph_div);
    this.sequence_div.css('top', top);

    // set widths of the central protein view
    var peptide_width = this.column3_div.outerWidth(true);
    var width = window_width - protein_list_width - peptide_width;
    set_outer_width(this.column2_div, width);
    width = get_content_width(this.column2_div);
    this.peptograph_canvas.set_width(width);
    set_outer_width(this.sequence_div, width);

    // move and resize the width right column 
    this.column3_div.css('left', window_width - peptide_width);

    // set the right column heights
    set_outer_height(this.column3_div, main_height);
    var height = this.column3_div.innerHeight();
    set_outer_height(this.peptide_list_div, Math.round(0.3*height));
    this.spectrum_canvas.set_height(Math.round(0.2*height));
    var top = get_bottom(this.peptide_info_div) + this.spectrum_div.outerHeight(true);
    set_outer_height(this.ion_table_div, height - top);
  }

  this.update = function() {
    this.resize_display();
    this.protein_list.highlight_selected_protein();
    this.peptograph_canvas.draw();
    this.spectrum_canvas.draw();
    this.sequence_view.update_panel();
    build_protein_info_panel(this.data, this.protein_info_div);  
    build_peptides_panel(this.data, this.peptide_list_div);
    build_peptide_info(this.data, this.peptide_info_div);
    build_ion_table_div(this.data, this.ion_table_div);
  }

  this.register_callbacks = function() {
    var _this = this;

    function onkeydown(event) {
      var c = String.fromCharCode(event.keyCode).toUpperCase();
      if (c == 'N') {
        var i = _this.protein_list.i_protein + 1;
        if (i < _this.data.sorted_seqids.length) {
          var seqid = _this.data.sorted_seqids[i];
          _this.data.controller.pick_protein_callback(seqid)();
          return false;
        }
      } else if (c == 'P') {
        var i = _this.protein_list.i_protein - 1;
        if (i >= 0) {
          var seqid = _this.data.sorted_seqids[i];
          _this.data.controller.pick_protein_callback(seqid)();
          return false;
        }
      } else if (c == 'J') {
        var protein = get_current_protein(_this.data);
        var i_source = protein.i_source_selected;
        var i_peptide = protein.i_peptide_selected + 1;
        var peptides = protein.sources[i_source].peptides;
        if (i_peptide < peptides.length) {
          _this.data.controller.pick_peptide_callback(protein, i_source, i_peptide)();
        }
      } else if (c == 'K') {
        var protein = get_current_protein(_this.data);
        var i_source = protein.i_source_selected;
        var i_peptide = protein.i_peptide_selected - 1;
        if (i_peptide >= 0) {
          _this.data.controller.pick_peptide_callback(protein, i_source, i_peptide)();
        }
      } else if (c == 'D') {
        var protein = get_current_protein(_this.data);
        var i_source = protein.i_source_view;
        if (i_source < protein.sources.length-1) {
          _this.data.controller.pick_slice(i_source+1);
          _this.data.observer();
          return false;
        }
      } else if (c == 'U') {
        var protein = get_current_protein(_this.data);
        var i_source = protein.i_source_view;
        if (i_source > 0) {
          _this.data.controller.pick_slice(i_source-1);
          _this.data.observer();
          return false;
        }
      };
      return;
    }

    this.data.observer = function() { _this.update(); };
    document.onkeydown = onkeydown;
    $(window).resize(_this.data.observer);
    this.data.observer();
  }

  this.init_data(data);
  this.init_widgets();
  this.register_callbacks();
}


$(function () { var app = new Pepto(data) });




