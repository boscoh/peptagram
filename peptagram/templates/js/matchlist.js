function count_matches_in_source(source) {
  source.attr = {}
  source.attr.n_match = 0;
  source.attr.n_match_unique = 0;
  var matches = source.matches;
  var n_match_in_slice = 0;
  var sequences = [];
  for (var i=0; i<matches.length; i++) {
    var match = matches[i];
    if (matches[i].attr.is_unique) {
      source.attr.n_match_unique += 1;
    }
    source.attr.n_match += 1;
    var seq = match.sequence;
    if (sequences.indexOf(seq) < 0) {
      sequences.push(seq);
    }
    source.attr.n_peptide = sequences.length;
  }
}


function build_matches_panel(data, div) {
  div.empty();
  var protein = data.controller.get_current_protein();
  var i_source = protein.i_source_selected;
  var source = protein.sources[i_source];
  count_matches_in_source(source);

  div.append('n_match: ' + source.attr.n_match + '<br>');
  div.append('n_peptide: ' + source.attr.n_peptide + '<br>');
  div.append('<br>');

  var matches = source.matches;

  if (matches.length == 0) {
    return;
  }
  
  var table = $('<table>');
  table.css('text-align', 'left');
  div.append(table);

  // heading
  var tr = $('<tr>')
  table.append(tr);
  tr.append($('<td>').text('#'));
  tr.append($('<td>').text('i'));
  tr.append($('<td>').text('seq'));

  for (i_match=0; i_match<matches.length; i_match++) {
    var match = matches[i_match];

    var tr = $('<tr>')
    table.append(tr);

    var td = $('<td>');
    td.append(i_match+1);
    tr.append(td);

    var td = $('<td>');
    td.append(match.i+1);
    tr.append(td);

    var td = $('<td>');
    var peptide_link = $("<a>");
    if ((i_match == protein.i_match_selected)) {
      peptide_link.addClass('highlight_peptide');
    }
    peptide_link.attr('href', '#');
    peptide_link.addClass('link');
    var seq = match['sequence'];
    var modified = [];
    var is_nterminal = false;
    var is_cterminal = false;
    if ('modifications' in match) {
      var modifications = match.modifications;
      for (var i=0; i<seq.length; i++) {
        modified.push(false);
      }
      for (var i=0; i<modifications.length; i++) {
        modified[modifications[i].i] = true;
        if (modifications[i].i == -1) {
          is_nterminal = true;
        }
        else if (modifications[i].i == seq.length) {
          is_cterminal = true;
        }
      }
    }
    if (is_nterminal) {
      peptide_link.append($('<span>').append('>').addClass('modification'));
    }
    for (var i=0; i<seq.length; i++) {
      if (modified[i]) {
        var s = $('<span>').append(seq[i]).addClass('modification');
      } else {
        var s = seq[i];    
      }
      peptide_link.append(s);
    }
    if (is_cterminal) {
      peptide_link.append($('<span>').append('<').addClass('modification'));
    }
    peptide_link.click(data.controller.pick_match_callback(protein, i_source, i_match));
    td.append(peptide_link);
    tr.append(td);
  }
}
