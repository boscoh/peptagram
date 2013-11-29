function build_peptides_panel(data, div) {
  div.empty();
  var protein = data.controller.get_current_protein();
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
    var modified = [];
    var is_nterminal = false;
    var is_cterminal = false;
    if ('modifications' in peptide.attr) {
      var modifications = peptide.attr.modifications;
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
    peptide_link.click(data.controller.pick_peptide_callback(protein, i_source, i_peptide));
    td.append(peptide_link);
    tr.append(td);
  }
}
