


function attr_empty(dict, key) {
  if (!(key in dict)) {
    return true;
  }
  if (dict[key].length == 0) {
    return true;
  }
  return false;
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
          val += ' '
        }
        var res_num = modifications[j].i + 1;
        val += res_num + ":" + modifications[j].mass;
      }  
    }
    s += key + ": " + val + "<br>";
  }
  return s;
}


function draw_highlight_box(canvas, x, y, w, h, color) {
  canvas.box(x-1, y-1, w+2, h+2, "white", 1);
  canvas.box(x-3, y-3, w+6, h+6, color, 3);
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


function set_max_outer_height(div, max_height) {
  var scroll_height = div[0].scrollHeight;
  var height = max_height;
  console.log(scroll_height, max_height);
  if (scroll_height < max_height) {
    height = scroll_height;
  }
  set_outer_height(div, height);
}


function set_max_outer_width(div, max_width) {
  var scroll_width = div[0].scrollWidth;
  var width = max_width;
  if (scroll_width < max_width) {
    width = scroll_width;
  }
  set_outer_width(div, width);
}


function get_outer_width(div) {
  return div.outerWidth(true);
}


function get_outer_height(div) {
  return div.outerHeight(true);
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


function get_right(div) {
  return div.position().left + div.outerWidth(true);
}


function set_top(div, top) {
  div.css('top', top);
}


function set_left(div, left) {
  div.css('left', left);
}

