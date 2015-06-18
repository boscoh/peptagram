
//
// peptagram.js a javascript app for displaying multiple 
// Label-Free protein identification using MS/MS analysis
// (c) 2013 Bosco Ho



function build_protein_info_panel(data, div) {
  div.empty();
  var protein = data.controller.get_current_protein();
  div.append(protein['description']);
  div.append("<br>----<br>");
  div.append(dict_html(protein['attr']));
}



// Converts intensities [-1, 1] into an RGB color string
function ColorIntensityPalette() {
  this.positive = [255, 0, 0];
  this.neutral = [200, 200, 200];
  this.negative = [0, 93, 144];

  this.rgb = function(color) {
      rgb = "rgb(" + 
            color[0] + "," +
            color[1] + "," +
            color[2] + ")";
      return rgb;
  }

  this.str = function(intensity) {
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
    var height = this.height;
    var half_height = height/2;
    var y_min = this.y;
    var x2 = this.x + this.bar_width;
    var y_max = height;
    for (var i=0; i<height; i++) {
      f = (half_height - i)/half_height;
      y = y_min + i;
      this.canvas.line(this.x, y, x2, y, this.palette.str(f), 1);
      if ((this.names.length == 2) && (i==half_height)) {
        break;
      }
    }
    var canvas = this.canvas;
    var write_name = function(txt, y) {
      canvas.text(txt, x2 + 2, y, this.font, "#999", "left");
    }
    write_name(this.names[0], this.y);
    write_name(this.names[1], this.y + 0.5*this.height);
    if (this.names.length>2) {
      write_name(this.names[2], this.y + this.height);
    }
  }
}


function PeptagramWidget(canvas, data, color_bar) {
  this.canvas = canvas;
  this.seqid = null;
  this.data = data;
  this.color_bar = color_bar;
  this.pressed = false;
  this.down = false;
  this.up = false; 
  this.x = 0;
  this.y = 0;
  this.protein = this.data.controller.get_current_protein();

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
    this.protein = this.data.controller.get_current_protein();

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

    // fit around the color_bar which is implicitly to the left
    this.x = this.color_bar.x + this.color_bar.width + 10 + padding;
    this.y = padding;
    this.draw_width = this.canvas.canvas_dom.width - this.x - name_width - freq_width;
    this.draw_height = this.canvas.canvas_dom.height - this.y - padding;
    var slice_height = this.get_diff_height(1, 0);
    var freq_height = slice_height/2;

    // draw background area
    this.canvas.solid_box(
        this.x, this.y, this.draw_width, this.draw_height, "white");

    // draw slice freq sidebar
    var max_n_pep = 0;
    for (var j=0; j<sources.length; j++) {
      if (sources[j].matches.length > max_n_pep) {
        max_n_pep = sources[j].matches.length;
      }
    }
    for (var i_source=0; i_source<sources.length; i_source++) {
      var n_pep = sources[i_source].matches.length;
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
        this.y_from_j(this.protein.i_source_selected),
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

    // draw matches
    for (var j=0; j<sources.length; j++) {
      var matches = sources[j].matches;
      for (var i=0; i<matches.length; i++) {
        var match = matches[i];
        if (match.intensity === "") {
          var color = "lightyellow";
        } else {
          var color = this.color_bar.palette.str(match.intensity);
        }
        this.canvas.solid_box(
            this.x_from_i(match.i), 
            this.y_from_j(j), 
            this.get_diff_width(match.i, match.j), 
            slice_height, 
            color);
      }
    }

    //  highlight sequence_view area
    draw_highlight_box(
        this.canvas, 
        this.x_from_i(this.protein.i_res_view), 
        this.y_from_j(this.protein.i_source_selected), 
        this.get_diff_width(0, this.data.n_res_in_view), 
        slice_height,
        "#AAA");

    //  highlight selected match
    var i_source_selected = this.protein.i_source_selected;
    var matches = sources[i_source_selected].matches;
    if (matches.length > 0) {
      if (this.protein.i_match_selected >= 0) {
        match = matches[this.protein.i_match_selected];
        draw_highlight_box(
            this.canvas, 
            this.x_from_i(match.i),
            this.y_from_j(i_source_selected),
            this.get_diff_width(match.i, match.j),
            slice_height,
            "#3F3");
      }
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

    // if click over match: select it
    var source = this.protein.sources[i_source];
    for (var i=0; i<source.matches.length; i++) {
      var match = source.matches[i];
      if ((match.i<=i_res) && (i_res<match.j)) {
        this.data.controller.pick_match(this.protein, i_source, i);
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

  this.sequence_number_line = function(n_res) {
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

  this.build = function() {
    this.seqid = this.data.selected_seqid;
    this.protein = this.data.controller.get_current_protein();
    this.i_source = this.protein.i_source_selected;
    this.i_match_selected = this.protein.i_match_selected;
   
    var matches = this.protein.sources[this.i_source].matches;
    var sequence = this.protein.sequence;
    var n_res = sequence.length;

    var is_peptides_intersect = function(i1, i2) {
      var pep1 = matches[i1];
      var pep2 = matches[i2];
      var not_intersect = ((pep1.j <= pep2.i) || (pep2.j <= pep1.i));
      return !not_intersect;
    }

    this.div.empty();

    var sequence_header = this.sequence_number_line(n_res);
    var pre = $('<pre>')
    this.div.append(pre);

    pre.append(sequence_header);
    pre.append('<br>');

    i_res = 0;
    for (var i_match=0; i_match<matches.length; i_match++) {
      match = matches[i_match];
      if (i_res < match.i) {
        pre.append(sequence.slice(i_res, match.i));
        i_res = match.i;
      }
      var peptide_link = $("<a>");
      if (this.i_match_selected >= 0) {
        if (this.i_match_selected < matches.length) {
          if (is_peptides_intersect(i_match, this.i_match_selected)) {
            peptide_link.addClass('highlight_peptide');
          }
        }
      }
      peptide_link.attr('href', '#');
      peptide_link.addClass('link');
      peptide_link.text(sequence.slice(i_res, match.j));
      peptide_link.click(this.data.controller.pick_match_callback(this.protein, this.i_source, i_match));
      pre.append(peptide_link);
      i_res = match.j;
    }
    if (i_res < n_res) {
      pre.append(sequence.slice(i_res, n_res));
    }
    this.data.n_res_in_view = Math.ceil(this.div.width()/this.char_width());
  }

  this.update_panel = function() {
    if ((this.seqid !== this.data.selected_seqid) ||
        (this.i_source !== this.protein.i_source_selected) ||
        (this.i_match_selected != this.protein.i_match_selected)) {
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


function build_match_info(data, div) {
  div.empty();
  var match = data.controller.get_selected_match();
  if (match != null) {
    div.append(dict_html(match['attr']));
    div.append("<br>");
  }
}


function Peptagram(data) {

  this.init_divs = function() {
    this.header_div = $('#header');
    this.column1_div = $('#column1');
    this.protein_control_div = $("#protein_control");
    this.protein_list_div = $("#protein_list");
    this.column2_div = $('#column2');
    this.protein_info_div = $('#protein_info')
    this.peptagram_div = $("#peptagram")
    this.sequence_div = $("#sequence");
    this.column3_div = $("#column3");
    this.peptide_list_div = $("#peptide_list");
    this.peptide_info_div = $('#peptide_info')
    this.ion_table_div = $('#ion_table')
    this.spectrum_div = $("#spectrum");
    this.peptide_match_info_div = $("#peptide_match_info");
    block_bounce_except_for_touchscroll();
    var _this = this;
    $(window).resize(function() { _this.update(); });
    this.update();
  }

  this.init_widgets = function() {
    this.header_div.text(this.data.title);
    document.title = this.data.title;

    // Build interactive page objects
    this.protein_list = new ProteinList(
        this.protein_control_div, this.protein_list_div, this.data);

    this.peptagram_div.css('width', '100%');
    this.peptagram_canvas = new CanvasWidget(this.peptagram_div, this.data.bg_color, true);
    var width = this.peptagram_canvas.div.width();
    this.color_bar = new ColorBarWidget(
        this.peptagram_canvas, this.data['color_names'], this.data.canvas_font);
    this.color_bar.x = 0;
    this.color_bar.y = 5;
    var peptagram_widget = new PeptagramWidget(
        this.peptagram_canvas, this.data, this.color_bar);
    peptagram_widget.x = 30;
    this.peptagram_canvas.push(this.color_bar);
    this.peptagram_canvas.push(peptagram_widget);

    this.sequence_view = new SequenceView(this.sequence_div, this.data);
    this.spectrum_canvas = new CanvasWidget(this.spectrum_div, this.data.bg_color, true);

    this.spectrum_widget = new SpectrumWidget(this.spectrum_canvas, this.data);
    this.spectrum_canvas.push(this.spectrum_widget);
    this.ion_table = new IonTable(this.ion_table_div, this.data);
  }

  this.resize_display = function() {
    var window_width = $(window).width();
    var window_height = $(window).height();

    // set columns at the right height
    var header_height = get_bottom(this.header_div);
    set_top(this.column1_div, header_height);
    set_top(this.column2_div, header_height);
    set_top(this.column3_div, header_height);
    var main_height = window_height - header_height;

    // put the protein widget bars in the right place
    set_outer_height(this.column1_div, main_height);
    var column1_height = get_content_height(this.column1_div);
    var top =  get_bottom(this.protein_control_div);
    set_outer_height(this.protein_list_div, column1_height - top);

    // resize column1_width if necessary
    var column1_width = get_outer_width(this.column1_div);
    if (window_width < 1600) {
      var default_width = 180;
    } else {
      var default_width = 450;
    }
    if (column1_width != default_width) {
      column1_width = default_width;
      set_outer_width(this.column1_div, column1_width);
    }

    // move the central column to its right place
    set_left(this.column2_div, column1_width);

    // set heights and tops of central protein view
    set_outer_height(this.column2_div, main_height);
    var top = get_bottom(this.protein_info_div);
    set_top(this.peptagram_div, top);

    // set the actual peptagram height
    var peptagram_height = get_content_height(this.column2_div) - 1;
    peptagram_height += - top - get_outer_height(this.sequence_div)
    if (this.short_peptagram_height < peptagram_height) {
      peptagram_height = this.short_peptagram_height;
    }
    set_outer_height(this.peptagram_div, peptagram_height);
    top = get_bottom(this.peptagram_div);
    set_top(this.sequence_div, top);

    // check column3 variable width
    var column3_width = get_outer_width(this.column3_div)
    if (window_width < 1600) {
      var default_width = 240;
    } else {
      var default_width = 450;
    }
    if (column3_width != default_width) {
      column3_width = default_width;
      set_outer_width(this.column3_div, column3_width);
    }

    // set widths of the central protein view
    var width = window_width - column1_width - column3_width;
    set_outer_width(this.column2_div, width);
    width = get_content_width(this.column2_div);
    set_outer_width(this.sequence_div, width);
    set_outer_width(this.peptagram_div, width);

    // move and resize the width right column 
    set_left(this.column3_div, window_width - column3_width);

    // set the right column heights
    set_outer_height(this.column3_div, main_height);
    var height = get_content_height(this.column3_div);
    set_outer_height(this.peptide_list_div, Math.round(0.2*height));
    var top = get_bottom(this.peptide_list_div);
    set_top(this.spectrum_div, top);
    set_outer_height(this.spectrum_div, Math.round(0.2*height));
    var top = get_bottom(this.peptide_info_div) + this.spectrum_div.outerHeight(true);
    set_outer_height(this.ion_table_div, height - top);
  }

  this.update = function() {
    if (this.data != null) {
      var n_source = this.data.controller.get_current_protein().sources.length;
      this.short_peptagram_height = 3000;
      if (n_source < 5) {
        this.short_peptagram_height = this.color_bar.height*n_source + 10;
      }
    }    

    this.resize_display();

    if (this.data != null) {
      this.protein_list.update();
      this.peptagram_canvas.draw();
      this.spectrum_widget.update();
      this.sequence_view.update_panel();
      build_protein_info_panel(this.data, this.protein_info_div);  
      build_matches_panel(this.data, this.peptide_list_div);
      build_match_info(this.data, this.peptide_info_div);
      this.ion_table.update();
    }
  }

  this.register_callbacks = function() {
    var _this = this;
    document.onkeydown = function(event) { 
      var c = String.fromCharCode(event.keyCode).toUpperCase();
      _this.data.controller.onkeydown(c); 
    }
    this.data.observer = function() { _this.update(); };
    this.data.observer();
  }

  this.init_data = function(data) {
    this.data = data;
    this.data.controller = new DataController(data);
    this.init_widgets();
    this.register_callbacks();
  }

  this.data = null;
  this.init_divs();
}


$(function () { 
  window.app = new Peptagram();
  window.load_data = function(data) { app.init_data(data); }
  load_script('data.jsonp');
});




