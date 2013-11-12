


// ProteinBarWidget draws a simple distribution of peptides in a
// protein using data.mask to mask certain peptides
function ProteinBarWidget(canvas, data, seqid) {
  this.canvas = canvas;
  this.seqid = seqid;
  this.data = data;
  this.protein = this.data.proteins[seqid];
  this.color = this.data.text_color; 
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
      var bg_color = this.data.bg_color;
    }
    this.canvas.solid_box(this.x, this.y, this.width, this.height, bg_color);
    this.canvas.solid_box(this.x, this.y + this.height/2, this.width, 1, '333');
    for (var j=0; j<this.protein.sources.length; j++) {
      var peptides = this.protein.sources[j].peptides;
      for (var i=0; i<peptides.length; i++) {
        var peptide = peptides[i];
        if (this.data.mask >= peptide.mask) {
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
    if (this.data.mask_labels.length > 0) {
      this.control_div.append('FPE: ');
      for (var i=0; i<this.data.mask_labels.length; i++) {
        var button = $("<input>");
        button.attr('name', 'mask');
        button.attr('type', 'radio');
        button.attr('value', this.data.mask_labels[i]);
        if (i == 0) {
          button.attr('checked', true);
        }
        this.control_div.append(button);
        this.control_div.append(this.data.mask_labels[i] + ' ');
      }
      this.control_div.append('<br>');
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

  this.update = function() {
    this.i_protein = this.data.controller.get_i_protein();
    if (this.i_old_protein == this.i_protein) {
      return;
    }
    this.data.selected_seqid = this.data.sorted_seqids[this.i_protein];
    if (this.i_old_protein != null) {
      this.protein_divs[this.i_old_protein].css('color', '#666');
      this.protein_widgets[this.i_old_protein].draw();
    }
    this.protein_divs[this.i_protein].css({'color':'#333'});
    this.protein_widgets[this.i_protein].draw();
    this.i_old_protein = this.i_protein;
  }


  this.build_list = function() {
    this.main_div.empty();
    while (this.protein_widgets.length > 0) {
      this.protein_widgets.pop();
    }
    while (this.protein_divs.length > 0) {
      this.protein_divs.pop();
    }

    this.main_div.append('<br>');
    this.build_sorting_msg();
    this.sort_key = this.attr_select.val();
    this.sort_direction = $("input[name=direction]:checked").attr('value');
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
        _this.update();
      }
    }
    progressive_build_list();
  }

  this.redraw_mask = function() {
    var val = $("input[name=mask]:checked").attr('value');
    this.data.mask = parseFloat(val);
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