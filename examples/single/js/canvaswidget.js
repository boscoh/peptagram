   ///////////////////////////////////////////////////
// Wrapper for the 2D canvas for HTML5
///////////////////////////////////////////////////


function is_canvas_suppported() {
  return ((!!window.HTMLCanvasElement) &&
          (!!window.CanvasRenderingContext2D));
}


// CanvasWidget is an object that interfaces between the 
// HTML5 canvas object and javascript. It wraps everything
// around a javascript object

// CanvasWidget contains a observer list: this.widgets.
// This is registered through this.push(widget)
// A widget has to implement:
//  - this.draw()
//  - this.pressed: boolean
//  - this.drag(x, y)
//  - this.x: Number
//  - this.y: Number
//  - this.width: Number
//  - this.height: Number
// The widget thus gets a simple IO with the user
// mainly through the this.drag function.

var CanvasWidget = function(div, bg_color) {

  this.pos_dom = function() {
    var curr_dom = this.canvas_dom;
    var curr_left = curr_top = 0;
    if (curr_dom.offsetParent) {
      curr_left = curr_dom.offsetLeft;
      curr_top = curr_dom.offsetTop;
      while (curr_dom = curr_dom.offsetParent) {
        curr_left += curr_dom.offsetLeft;
        curr_top += curr_dom.offsetTop;
      }
    }
    curr_dom = this.canvas_dom;
    do {
      curr_left -= curr_dom.scrollLeft || 0;
      curr_top -= curr_dom.scrollTop || 0;
    } while (curr_dom = curr_dom.parentNode);
    return [curr_left, curr_top];
  }

  this.set_width = function(w) { 
    this.canvas_dom.width = w;
    this.div.width(w);
  }

  this.set_height = function(h) {
    this.canvas_dom.height = h; 
    this.div.height(h);
  }

  this.x = function() { return this.pos_dom(this.canvas_dom)[0]; }

  this.y = function() { return this.pos_dom(this.canvas_dom)[1]; }

  this.width = function() { return this.canvas_dom.width; }

  this.height = function() { return this.canvas_dom.height; }

  this.set_scale = function() {
    var h = this.height();
    var w = this.width();
    this.scale = Math.sqrt(h*h + w*w);
  }

  this.extract_mouse_xy = function(event) {
    this.x_mouse = event.clientX - this.x();
    this.y_mouse = event.clientY - this.y();
    if (event.touches) {
      this.x_mouse = event.touches[0].clientX - x;
      this.y_mouse = event.touches[0].clientY - y;
    }
  }
  
  this.line = function(x1, y1, x2, y2, color, width) {
    this.draw_context.beginPath();
    this.draw_context.moveTo(x1, y1);
    this.draw_context.lineTo(x2, y2);
    this.draw_context.closePath();
    this.draw_context.lineWidth = width;
    this.draw_context.strokeStyle = color;
    this.draw_context.stroke();

  };

  this.box = function(x, y, w, h, color, width) {
    this.line(x, y, x + w, y, color, width);
    this.line(x + w, y, x + w, y + h, color, width);
    this.line(x + w, y + h, x, y + h, color, width);
    this.line(x, y + h, x, y, color, width);
  }


  this.solid_box = function(x, y, w, h, color) {
    this.draw_context.fillStyle = color;
    this.draw_context.fillRect(x, y, w, h);
  }

  this.solid_circle = function(x, y, r, color, edgecolor) {
    this.draw_context.beginPath();
    this.draw_context.arc(x, y, r, 0, 2*Math.PI, true);
    this.draw_context.closePath();
    this.draw_context.fillStyle = color;
    this.draw_context.fill();
  };

  this.circle = function(x, y, r, color, width) {
    this.draw_context.beginPath();
    this.draw_context.arc(x, y, r, 0, 2*Math.PI, true);
    this.draw_context.closePath();
    this.draw_context.strokeStyle = color;
    this.draw_context.lineWidth = width;
    this.draw_context.stroke();
  };

  this.solid_line = function(x1, y1, x2, y2, th, color, edgecolor) {
    var dx = y1 - y2;
    var dy = x2 - x1;
    var d  = Math.sqrt(dx*dx + dy*dy);
    dx /= d;
    dy /= d;
    this.draw_context.beginPath();
    this.draw_context.moveTo(x1 + dx * th, y1 + dy * th);
    this.draw_context.lineTo(x2 + dx * th, y2 + dy * th);
    this.draw_context.lineTo(x2 - dx * th, y2 - dy * th);
    this.draw_context.lineTo(x1 - dx * th, y1 - dy * th);
    this.draw_context.closePath();
    this.draw_context.fillStyle = color;
    this.draw_context.fill();
    this.draw_context.strokeStyle = edgecolor;
    this.draw_context.lineWidth = 1;
    this.draw_context.stroke();
  };

  this.quad = function(
      x1, y1, x2, y2, x3, y3, x4, y4, 
      color, edgecolor) {
    this.draw_context.beginPath();
    this.draw_context.moveTo(x1, y1);
    this.draw_context.lineTo(x2, y2);
    this.draw_context.lineTo(x3, y3);
    this.draw_context.lineTo(x4, y4);
    this.draw_context.closePath();
    this.draw_context.fillStyle = color;
    this.draw_context.fill();
    this.draw_context.strokeStyle = edgecolor;
    this.draw_context.lineWidth = 1;
    this.draw_context.stroke();
  };

  this.text = function(text, x, y, font, color, align) {
    this.draw_context.fillStyle = color;
    this.draw_context.lineStyle = color;
    this.draw_context.lineWidth = 0.5;
    this.draw_context.font = font;
    if (typeof align === 'undefined') {
      var align = 'left';
    }
    this.draw_context.textAlign = align;
    this.draw_context.textBaseline = 'middle';
    this.draw_context.fillText(text, x, y);
  }

  this.get_textwidth = function(text, font) {
    this.draw_context.font = font;
    this.draw_context.textAlign = 'center';
    return this.draw_context.measureText(text).width;
  }
  
  this.draw_popup = function(x, y, text, fillstyle) {
    var h = 20;
    var w = this.get_textwidth(text, '10px sans-serif') + 20;
    var y1 = 30;
    var arrow_w = 5;
    this.draw_context.beginPath();
    this.draw_context.moveTo(x-arrow_w, y-y1);
    this.draw_context.lineTo(x, y);
    this.draw_context.lineTo(x+arrow_w, y-y1);
    this.draw_context.lineTo(x+w/2, y-y1);
    this.draw_context.lineTo(x+w/2, y-h-y1);
    this.draw_context.lineTo(x-w/2, y-h-y1);
    this.draw_context.lineTo(x-w/2, y-y1);
    this.draw_context.closePath();
    this.draw_context.fillStyle = fillstyle;
    this.draw_context.fill();
    this.text(
        text, x, y-h/2-y1, '10px sans-serif', '#000', 'center');
  };

  this.draw_background = function() {
    this.width = this.canvas_dom.width;
    this.height = this.canvas_dom.height;
    this.half_width = this.width/2;
    this.half_height = this.height/2;
    var w = this.width;
    var h = this.height;
    this.scale = Math.sqrt(h*h + w*w);
    this.draw_context.clearRect(0, 0, w, h); 
  }


  this.widgets = [];

  this.push = function(widget) {
    this.widgets.push(widget);
  }

  this.draw = function() {
    this.draw_background();
    for (var i=0; i<this.widgets.length; i++) {
      this.widgets[i].draw();
    }
  }

  this.in_widget = function(widget, x, y) {
    if ((x < widget.x) || (x > (widget.x + widget.width))) {
      return false;
    }
    if ((y < widget.y) || (y > (widget.y + widget.height))) {
      return false;
    }
    return true;
  }

  this.mousedrag = function(x, y) {
    var something_happened = false;
    for (var i=0; i<this.widgets.length; i++) {
      var widget = this.widgets[i];
      if (this.in_widget(widget, x, y)) {
        if ('drag' in widget) {
          widget.drag(x, y);
          something_happened = true;
        }
      }
    }
    // if (something_happened) {
    //   this.draw();
    // }
  }

  this.mousedown = function(event) {
    this.extract_mouse_xy(event);
    var x = this.x_mouse;
    var y = this.y_mouse;
    for (var i=0; i<this.widgets.length; i++) {
      var widget = this.widgets[i];
     if (this.in_widget(widget, x, y)) {
        widget.pressed = true;
        widget.down = true;
        widget.up = false;
      }
    }
    this.mousedrag(x, y);
    for (var i=0; i<this.widgets.length; i++) {
      var widget = this.widgets[i];
     if (this.in_widget(widget, x, y)) {
        widget.down = false;
      }
    }
    event.preventDefault();
  }

  this.mousemove = function(event) {
    this.extract_mouse_xy(event);
    var x = this.x_mouse;
    var y = this.y_mouse;
    for (var i=0; i<this.widgets.length; i++) {
      var widget = this.widgets[i];
     if (this.in_widget(widget, x, y)) {
        widget.down = false;
        widget.up = false;
      }
    }
    this.mousedrag(x, y);
    event.preventDefault();
  }

  this.mouseup = function(event) {
    this.extract_mouse_xy(event);
    var x = this.x_mouse;
    var y = this.y_mouse;
    for (var i=0; i<this.widgets.length; i++) {
      var widget = this.widgets[i];
      if (this.in_widget(widget, x, y)) {
        widget.pressed = false;
        widget.down = false;
        widget.up = true;
      }
    }
    this.mousedrag(x, y);
    event.preventDefault();
  }

  this.div = div;
  if (!is_canvas_suppported()) {
    this.div.append('Sorry, no canvs2d detected in your browser!');
    return;
  }

  this.x_mouse;
  this.y_mouse;
  this.canvas = $('<canvas>').css('background-color', bg_color);
  this.div.append(this.canvas);
  this.canvas_dom = this.canvas[0];
  this.draw_context = this.canvas_dom.getContext('2d');
  this.set_width(this.div.width());
  this.set_height(this.div.height());
  this.set_scale();

  var _this = this;
  this.canvas_dom.addEventListener(
    'mousedown',
    function(e) { _this.mousedown(e); }, 
    false);
  this.canvas_dom.addEventListener(
    'mouseup',
    function(e) { _this.mouseup(e); }, 
    false);
  this.canvas_dom.addEventListener(
    'mousemove',
    function(e) { _this.mousemove(e); }, 
    false);
}






