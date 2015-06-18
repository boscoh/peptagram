


// ProteinBarWidget draws a simple distribution of matches
// in a protein
function ProteinBarWidget(canvas, data, seqid) {
  this.canvas = canvas;
  this.seqid = seqid;
  this.data = data;
  this.protein = this.data.proteins[seqid];
  this.color = this.data.text_color; 
  this.x = 0;
  this.y = 0;

  this.x_from_i = function(i) {
    return i/this.protein.length*this.width + this.x;
  }

  this.draw = function() {
    this.width = this.canvas.width;
    this.height = this.canvas.height;

    if (this.seqid == this.data.selected_seqid) {
      var bg_color = this.data.select_bg_color;
    } else {
      var bg_color = this.data.bg_color;
    }
    this.canvas.solid_box(this.x, this.y, this.width, this.height, bg_color);
    this.canvas.solid_box(this.x, this.y + this.height/2, this.width, 1, '#333');
    for (var j=0; j<this.protein.sources.length; j++) {
      var matches = this.protein.sources[j].matches;
      for (var i=0; i<matches.length; i++) {
        var match = matches[i];
        var x = this.x_from_i(match.i);
        var w = this.x_from_i(match.j) - x;
        this.canvas.solid_box(x, this.y, w, this.height, this.color);
      }
    }
  }

  this.drag = function(x, y) {
    if (this.canvas.touch) {
      if (this.selected_seqid != this.seqid) {
        this.data.controller.pick_protein(this.seqid);
        this.data.observer();
      }  
    }
  }
}


function count_matches(protein) {
  protein.attr.n_match = 0;
  protein.attr.n_match_unique = 0;
  protein.attr.n_slice_populated = 0;
  var sources = protein.sources;
  for (var j=0; j<sources.length; j++) {
    var matches = sources[j].matches;
    var n_match_in_slice = 0;
    for (var i=0; i<matches.length; i++) {
      var match = matches[i];
      if (matches[i].attr.is_unique) {
        protein.attr.n_match_unique += 1;
      }
      protein.attr.n_match += 1;
      n_match_in_slice += 1;
    }
    if (n_match_in_slice > 0) {
      protein.attr.n_slice_populated += 1;
    }
  }
}


function ProteinList(control_div, column1_div, data) {
  this.main_div = column1_div
  this.control_div = control_div;
  this.data = data;
  this.protein_canvases = [];
  this.protein_divs = [];
  this.sort_key;
  this.sorting_msg_div;
  this.attr_select;

  this.build_controls = function() {
    // selectors for sorting protein attributes
    this.control_div.append("Sort by:");
    var keys = sorted_keys(this.data.controller.get_current_protein().attr);
    this.attr_select = $("<select>");
    for (i=0; i<keys.length; i++) {
      var key = keys[i];
      var option = $('<option>');
      option.attr('value', key);
      option.text(key);
      if (key == 'n_match') {
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

  this.build_list = function() {
    this.main_div.empty();
    while (this.protein_canvases.length > 0) {
      this.protein_canvases.pop();
    }
    while (this.protein_divs.length > 0) {
      this.protein_divs.pop();
    }

    this.main_div.append('<br>');
    this.build_sorting_msg();
    this.sort_key = this.attr_select.val();
    this.sort_direction = $("input[name=direction]:checked").attr('value');

    for (var seqid in this.data.proteins) {
      var protein = this.data.proteins[seqid];
      count_matches(protein);
    }

    this.data.controller.calc_sorted_seqids(
        this.sort_key, this.sort_direction);

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
        var val = protein.attr[_this.sort_key];
        var description = "" + (i_seqid+1) + ":";
        description += "[" + val + "] ";
        description += seqid + ": ";
        description += protein.description;

        var protein_div = $("<div>");
        protein_div.append($("<div>").text(description));
        protein_div.click(_this.data.controller.pick_protein_callback(seqid));
        _this.protein_divs.push(protein_div);
        _this.main_div.append(protein_div);

        var canvas_div = $("<div>");
        canvas_div.addClass("protein_bar");
        canvas_div.css('width', _this.main_div.width() - 10);
        canvas_div.css('height', 12);
        canvas_div.css('margin', '5px 0');
        canvas_div.css('width','100%');
        protein_div.append(canvas_div);

        var canvas = new CanvasWidget(canvas_div, "#EFEFEF");
        var protein_bar_widget = new ProteinBarWidget(canvas, data, seqid);
        _this.protein_canvases.push(canvas);

        canvas.push(protein_bar_widget);
        canvas.draw();
      }
      
      if (j_seqid < _this.data.sorted_seqids.length) {
        setTimeout(progressive_build_list, ms_wait);
      } else {
        _this.main_div.append('<br><br>')
        _this.sorting_msg_div.remove();
        _this.update();
      }
    }
    progressive_build_list();
  }

  this.redraw_bars = function() {
    var width = get_content_width(this.control_div);
    set_outer_width(this.sorting_msg_div, width);
    for (var i=0; i<this.protein_canvases.length; i++) {
      this.protein_canvases[i].draw();
    }
  }

  this.highlight_selected = function() {
    this.i_protein = this.data.controller.get_i_protein();
    if (this.i_old_protein == this.i_protein) {
      return;
    }
    this.data.selected_seqid = this.data.sorted_seqids[this.i_protein];
    if (this.i_old_protein != null) {
      if (this.i_old_protein >= this.protein_divs.length) {
        return;
      }
      this.protein_divs[this.i_old_protein].css('color', '#666');
      this.protein_canvases[this.i_old_protein].draw();
    }
      if (this.i_protein >= this.protein_divs.length) {
        return;
      }
    this.protein_divs[this.i_protein].css({'color':'#333'});
    this.protein_canvases[this.i_protein].draw();
    this.i_old_protein = this.i_protein;
  }

  this.update = function() {
    this.redraw_bars();
    this.highlight_selected();
  }

  this.register_callbacks = function() {
    var _this = this;
    this.attr_select.change(
        function() { _this.build_list(); });
    $("input[name=direction]").change(
        function() { _this.build_list(); });
  }

  this.build_controls();
  this.build_list();
  this.register_callbacks();
}