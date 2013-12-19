
//
// Pepto.js a javascript app for displaying multiple 
// Label-Free protein identification using MS/MS analysis
// (c) 2013 Bosco Ho


function build_protein_info_panel(data, div) {
  div.empty();
  var protein = data.controller.get_current_protein();
  div.append(protein['description']);
  div.append("<br>----<br>");
  div.append(dict_html(protein['attr']));
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

  this.is_peptides_intersect = function(i1, i2) {
    var pep1 = this.peptides[i1];
    var pep2 = this.peptides[i2];
    var not_intersect = ((pep1.j <= pep2.i) || (pep2.j <= pep1.i));
    return !not_intersect;
  }

  this.get_i_peptide = function(i_res) {
    for (var i_peptide=0; i_peptide<this.peptides.length; i_peptide++) {
      if ((this.peptides[i_peptide].i <= i_res) && 
          (i_res < this.peptides[i_peptide].j)) {
        return i_peptide;
      }
    }
    return -1;
  }

  this.make_link = function(i_peptide) {
    var peptide_link = $("<a>");
    if (this.i_peptide_selected < this.peptides.length) {
      if (this.is_peptides_intersect(i_peptide, this.i_peptide_selected)) {
        peptide_link.addClass('highlight_peptide');
      }
    }
    peptide_link.attr('href', '#');
    peptide_link.addClass('sequence_link');
    peptide_link.click(this.data.controller.pick_peptide_callback(this.protein, this.i_source, i_peptide));
    return peptide_link;
  }
  
  this.build = function() {
    this.seqid = this.data.selected_seqid;
    this.protein = this.data.controller.get_current_protein();
    this.i_source = this.protein.i_source_view;
    this.i_peptide_selected = this.protein.i_peptide_selected;
   
    this.div.empty();

    var pre = $('<pre>')
    this.div.append(pre);

    this.peptides = this.protein.sources[this.i_source].peptides;
    var sequence = this.protein.sequence;
    var n_res = sequence.length;

    var target = pre;
    current_i_peptide = -1;
    for (var i_res=0; i_res<sequence.length; i_res++) {
      // if (this.data.mask > peptide.mask) {
      //   continue;
      // }
      var i_peptide = this.get_i_peptide(i_res);
      if ((i_peptide >= 0) && (i_peptide != current_i_peptide)) {
        var peptide = this.peptides[i_peptide];
        var peptide_link = this.make_link(i_peptide);
        pre.append(peptide_link);
        target = peptide_link;
        current_i_peptide = i_peptide;
      }
      if ((i_peptide == -1) && (current_i_peptide >= 0)) {
        current_i_peptide = -1;
        target = pre;
      }
      if (i_res % 50 == 0) {
        if (i_res > 0) {
          pre.append('&nbsp; &nbsp;&nbsp;<br>');
        }
        num_str = '' + i_res;
        while (num_str.length < 5) {
          num_str = ' ' + num_str;
        }
        pre.append(num_str + ' ');
        if (current_i_peptide >= 0) {
          var peptide_link = this.make_link(i_peptide);
          pre.append(peptide_link);
          target = peptide_link;
        }
      } else if (i_res % 10 == 0) {
        target.append(' ');
      }
      target.append(sequence[i_res]);
    }
    this.data.n_res_in_view = Math.ceil(this.div.width()/this.char_width());
  }

  this.update = function() {
    if ((this.seqid !== this.data.selected_seqid) ||
        (this.i_source !== this.protein.i_source_view) ||
        (this.i_peptide_selected != this.protein.i_peptide_selected)) {
      this.build();
    }
  }

}


function build_peptide_info_panel(data, div) {
  div.empty();
  var peptide = data.controller.get_selected_peptide();
  div.append(dict_html(peptide.attr));
  div.append("<br>");
}



function Pepto(data) {
  this.init_data = function(data) {
    this.data = data;
    this.data.controller = new DataController(this.data);
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
    this.peptide_list_title_div = $("#peptide_list_title");
    this.peptide_list_div = $("#peptide_list");
    this.peptide_info_div = $('#peptide_info')
    this.ion_table_div = $('#ion_table')
    this.spectrum_div = $("#spectrum");
    this.peptide_match_info_div = $("#peptide_match_info");

    this.header_div.html(this.data.title);

    // Build interactive page objects
    this.protein_list = new ProteinList(
        this.protein_control_div, this.protein_list_div, this.data);
    this.sequence_view = new SequenceView(this.sequence_div, this.data);
    this.spectrum_canvas = new CanvasWidget(this.spectrum_div, "#EFEFEF");
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

    // move the central column to its right place
    var column1_width = get_outer_width(this.column1_div);
    set_left(this.column2_div, column1_width);
    var width = get_outer_width(this.sequence_div);
    this.column2_div.width(width);
    var width = get_content_width(this.column2_div);
    set_outer_width(this.peptide_list_div, width);
    set_outer_width(this.peptide_list_title_div, width);

    // set heights and tops of central column
    set_outer_height(this.column2_div, main_height);
    var column2_height = get_content_height(this.column2_div);
    var top = get_bottom(this.protein_info_div);
    set_top(this.sequence_div, top);
    set_outer_height(this.sequence_div, Math.round(0.4*column2_height));
    var top = get_bottom(this.sequence_div);
    set_top(this.peptide_list_title_div, top);
    var top = get_bottom(this.peptide_list_title_div);
    set_top(this.peptide_list_div, top);
    set_outer_height(this.peptide_list_div, column2_height - top);

    // move and resize the width right column 
    var left = get_right(this.column2_div);
    set_left(this.column3_div, left);
    var width = window_width - left;
    set_outer_width(this.column3_div, width);
    this.spectrum_canvas.set_width(width - 20);

    // set the right column heights
    set_outer_height(this.column3_div, main_height);
    var height = get_content_height(this.column3_div);
    this.spectrum_canvas.set_height(Math.round(0.3*height));
    var top = get_bottom(this.peptide_info_div) + this.spectrum_div.outerHeight(true);
    set_outer_height(this.ion_table_div, height - top);
  }

  this.update = function() {
    build_protein_info_panel(this.data, this.protein_info_div);  
    build_peptides_panel(this.data, this.peptide_list_div);
    build_peptide_info_panel(this.data, this.peptide_info_div);
    this.sequence_view.update();
    this.resize_display();
    // canvas drawing done after resize_display!
    this.protein_list.update();
    this.spectrum_widget.update();
    this.ion_table.update();
  }

  this.register_callbacks = function() {
    var _this = this;
    this.data.observer = function() { _this.update(); };
    document.onkeydown = function(event) { 
      var c = String.fromCharCode(event.keyCode).toUpperCase();
      _this.data.controller.onkeydown(c); 
    }
    $(window).resize(_this.data.observer);
    this.data.observer();
  }

  this.init_data(data);
  this.init_widgets();
  this.register_callbacks();
}


$(function () { var app = new Pepto(data) });




