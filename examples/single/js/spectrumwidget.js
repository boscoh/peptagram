function SpectrumWidget(canvas, data) {
  this.canvas = canvas;
  this.data = data;
  this.x = 0;
  this.y = 0;
  this.offset = 30;
  this.pressed = false;
  this.down = false;
  this.up = false; 
  this.mouse_in = false;
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

  this.update = function() {
    this.width = this.canvas.canvas_dom.width;
    this.height = this.canvas.canvas_dom.height;
    this.draw_width = this.width - 2*this.offset;
    this.draw_height = this.height - 2*this.offset;
    var label_color = '#BCB';

    var peptide = this.data.controller.get_selected_peptide();
    if (attr_empty(peptide, 'spectrum')) {
      this.canvas.div.css('display', 'none');
      return;
    } 
    this.canvas.div.css('display', 'block');

    this.spectrum = this.data.controller.get_labeled_spectrum();
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
      this.data.canvas_font, label_color, 'center')

    // draw zoom help line
    if (!this.data.zoom && this.mouse_in) {
      var i = this.data.i_peak;
      if (i) {
        var x1 = this.x_from_m(this.spectrum[i][0]);
         if ((x1 >= this.x + this.offset) &&
            (x1 <= this.x + this.offset + this.draw_width)) {
          var y1 = this.y + this.offset;
          var y2 = this.y + this.offset + this.draw_height;
          canvas.line(x1, y1, x1, y2, label_color, 1);
        }
      }
    }

    // draw peaks and peak labels
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
        color = "#DDD";
      }
      if (this.data.zoom) {
        canvas.line(x1, y1, x1, y2, color, 5);
      } else {
        canvas.line(x1, y1, x1, y2, color, 3);
      }
      if (label) {
        canvas.text(label, x1, y1 - 8, this.data.canvas_font, color);
      }
    }

    // draw frame
    canvas.line(
        this.x + this.offset,
        this.y + this.offset,
        this.x + this.offset,
        this.y + this.offset + this.draw_height,
        label_color, 1);
    canvas.line(
        this.x + this.offset,
        this.y + this.offset + this.draw_height,
        this.x + this.offset + this.draw_width,
        this.y + this.offset + this.draw_height,
        label_color, 1);

    // draw min, max values of mass/charage
    canvas.text(
      this.min_m.toFixed(0), 
      this.x + this.offset, 
      this.y + this.offset + this.draw_height + 10, 
      this.data.canvas_font, label_color, 'center')
    canvas.text(
      this.max_m.toFixed(0), 
      this.x + this.offset + this.draw_width, 
      this.y + this.offset + this.draw_height + 10, 
      this.data.canvas_font, label_color, 'center')
    if (this.data.zoom) {
      var label = 'm/z=' + peak[0].toFixed(0);
    } else {
      var label = 'm/z';
    }
    canvas.text(
      label, 
      this.x + this.offset + this.draw_width/2, 
      this.y + this.offset + this.draw_height + 10, 
      this.data.canvas_font, label_color, 'center')

    // draw max intensity label
    canvas.text(
      this.max_i.toFixed(0), 
      this.x + 2,
      this.y + this.offset - 8,
      this.data.canvas_font, label_color, 'left')
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
    var peptide = this.data.controller.get_selected_peptide();
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

  var _this = this;
  this.canvas.canvas.mouseout(function() { 
      _this.mouse_in = false;
      _this.data.observer();
  })

  this.canvas.canvas.mouseenter(function() { 
      _this.mouse_in = true;
      _this.data.observer();
  })
}



function IonTable(div, data) {
  this.data = data;
  this.div = div;

  this.make_td = function(s) {
    return $("<td>").append(s);
  }

  this.make_ion_heading = function(ion_type) {
    var s = ion_type;
    if (s.length == 1) {
      s = '&nbsp;' + s + '&nbsp;';
    }
    var a = $('<a>').html(s).addClass('link').attr('href','')
    var td = this.make_td(a);
    a.click(this.data.controller.toggle_ion_callback(ion_type));
    if (this.data.ion_types[ion_type]) {
      a.addClass('highlight_peptide');
    }
    return td;
  }

  this.mass_td = function(num, charge, ion_type) {
    var label = ion_label(num, charge, ion_type);
    var s = '';
    for (var i=0; i<this.spectrum.length; i++) {
      if (this.spectrum[i][2] === label) {
        s = $('<a>').addClass('link');
        s.text(this.spectrum[i][this.i_mass].toFixed(0));
        s.attr('href', '');
        s.click(this.data.controller.pick_peak_callback(i))
        break;
      }
    }
    return this.make_td(s);
  }

  this.change_mass_unit = function() {
    var mass_unit = $("input[name=mass_unit]:checked").attr('value');
    this.data.mass_unit = mass_unit;
    this.data.observer(); 
  }


  this.update = function() {
    this.div.empty();
    this.div.append('<br>')


    var peptide = this.data.controller.get_selected_peptide();
    if (attr_empty(peptide, 'spectrum')) {
      return;
    }

    for (var i=0, mass_unit; mass_unit=this.data.mass_units[i]; i++) {
      var button = $("<input>");
      button.attr('name', 'mass_unit');
      button.attr('type', 'radio');
      button.attr('value', mass_unit);
      if (data.mass_unit == mass_unit) {
        button.attr('checked', true);
      }
      div.append(button)
      div.append(mass_unit + '&nbsp;');
    }

    var _this = this;
    $("input[name=mass_unit]").change(function() { _this.change_mass_unit(); }) 

    var mass_unit = $("input[name=mass_unit]:checked").attr('value');
    if (mass_unit == this.data.mass_units[0]) {
      this.i_mass = 0;
    } else {
      this.i_mass = 4
    }

    var table = $("<table>");

    // build table header
    var tr = $("<tr>");
    tr.append(this.make_ion_heading("b(3+)"));
    tr.append(this.make_ion_heading("b(2+)"));
    tr.append(this.make_ion_heading("b"));
    tr.append(this.make_td("C"));
    tr.append(this.make_ion_heading("y"));
    tr.append(this.make_ion_heading("y(2+)"));
    tr.append(this.make_ion_heading("y(3+)"));
    table.append(tr);

    // build ion table
    this.spectrum = this.data.controller.get_labeled_spectrum();
    for (var i=0; i<peptide.sequence.length; i++) {
      var tr = $("<tr>");
      var num = i + 1;
      tr.append(this.mass_td(num, 3, "b"));
      tr.append(this.mass_td(num, 2, "b"));
      tr.append(this.mass_td(num, 1, "b"));
      pre = '' + num;
      if (pre.length == 1) {
        pre = '.' + pre;
      }
      post = '' + (peptide.sequence.length - i);
      if (post.length == 1) {
        post = post + '.';
      }
      tr.append(this.make_td(pre + '.' + peptide.sequence[i] + '.' + post));
      var num = peptide.sequence.length - i;
      tr.append(this.mass_td(num, 1, "y"));
      tr.append(this.mass_td(num, 2, "y"));
      tr.append(this.mass_td(num, 3, "y"));
      table.append(tr);
    }
    div.append(table);
    table.addClass('dashed_border');
  }

}


