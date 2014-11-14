
import sys, copy, os, re, datetime

from collections import OrderedDict


try:
    from reportlab import rl_config
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.platypus import BaseDocTemplate, Frame, PageTemplate, Spacer, Paragraph, NextPageTemplate, PageBreak, Table, TableStyle
    from reportlab.lib.styles import ParagraphStyle  
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.graphics.shapes import colors, Drawing, Rect, Line, Group, String
    from reportlab.lib.units import inch
    from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
    from reportlab.platypus.tableofcontents import TableOfContents  
except:
    raise ImportError("Couldn't import reportlab")

import peptidemass



# suppress creation timestamp and randomized document ID
rl_config.invariant = 1



# define styles
title_style = ParagraphStyle(
    name='Title', fontSize=24, leading=24, textColor=colors.Color(0.8, 0., 0.),
    alignment=TA_LEFT, fontName='Helvetica-Bold')
subtitle_style = ParagraphStyle(
    name='Subtitle', fontSize=12, leading=12, textColor='black',
    alignment=TA_LEFT, fontName='Helvetica-Bold')
bullet_style = ParagraphStyle(
    name='Bullet', fontSize=12, leading=12, spaceAfter=6, textColor='black',
    alignment=TA_LEFT, bulletIndent=0, leftIndent=10)
para_style = ParagraphStyle(
    name='Paragraph', leftIndent=30, spaceBefore=0, 
    spaceAfter=0, alignment=TA_JUSTIFY, fontSize=6, 
    leading=6)
heading1_style = ParagraphStyle(
    name='Heading1', fontSize=10, leading=6, 
    alignment=TA_LEFT, leftIndent=30, fontName='Helvetica-Bold')



def blank_page_draw(canvas, doc):
    canvas.saveState()
    canvas.restoreState()



def numbered_page_draw(canvas, doc):
    font_size = 6
    footer_bottom = doc.bottomMargin
    top_line = footer_bottom + font_size
    line_length = doc.width + doc.leftMargin

    canvas.saveState()
    canvas.setFont('Helvetica', font_size)
    text = "Page %d." % doc.page
    canvas.drawString(inch, footer_bottom, text)

    text = "made with http://boscoh.github.io/peptagram."
    canvas.drawRightString(line_length, footer_bottom, text)
    rect_left = line_length - canvas.stringWidth(text)
    link_rect = (line_length, footer_bottom, rect_left, top_line)
    link = "http://boscoh.github.io/peptagram"
    canvas.linkURL(link, link_rect)

    canvas.restoreState()


class TocDocTemplate(BaseDocTemplate):
    """
    A Doc Template that creates a TOC entry for every
    paragraph with style 'Heading1' and 'Heading2'
    """

    def afterFlowable(self, flowable):  
        "Registers TOC entries."  
        if flowable.__class__.__name__ == 'Paragraph':  
            text = flowable.getPlainText()  
            style = flowable.style.name  
            if style == 'Heading1':  
                level = 1  
            elif style == 'Heading2':  
                level = 2  
            else:  
                return  
            E = [level, text, self.page]  
            #if we have a bookmark name append that to our notify data  
            bn = getattr(flowable,'_bookmarkName',None)  
            if bn is not None: E.append(bn)  
            self.notify('TOCEntry', tuple(E))  


  
class ReportLabDoc():
    def __init__(self, pdf):
        self.elements = []
        self.chapter_num = 1
        self.doc = TocDocTemplate(pdf)
        self.page_templates = OrderedDict()
        self.build_page_templates()
        self.toc = TableOfContents()
        self.labels = ['b(3+)', 'b(2+)', 'b', '>', 'aa', '<', 'y', 'y(2+)', 'y(3+)']

    def add_toc(self):
        self.elements.append(self.toc)

    def build_page_templates(self):
        # normal frame as for SimpleFlowDocument
        frameT = Frame(
            self.doc.leftMargin, self.doc.bottomMargin, self.doc.width, 
            self.doc.height, id='normal')
        self.page_templates['First'] = \
            PageTemplate(
                id='First', frames=frameT, onPage=blank_page_draw)
        self.page_templates['OneCol'] = \
            PageTemplate(
                id='OneCol', frames=frameT, onPage=numbered_page_draw)

        # Two Columns
        frame1 = Frame(
            self.doc.leftMargin, self.doc.bottomMargin, 
            self.doc.width/2-6,  self.doc.height, 
            id='col1')
        frame2 = Frame(
            self.doc.leftMargin+self.doc.width/2+6, 
            self.doc.bottomMargin, 
            self.doc.width/2-6, 
            self.doc.height, 
            id='col2')
        self.page_templates['TwoCol'] = \
            PageTemplate(
                id='TwoCol',frames=[frame1,frame2], onPage=numbered_page_draw)

        self.doc.addPageTemplates(self.page_templates.values())
        
    def add_page_break(self):
        self.elements.append(PageBreak())

    def add_paragraph(self, txt, style=para_style):
        self.elements.append(Paragraph(txt, style))

    def add_bullet(self, txt, style=bullet_style):
        self.elements.append(
            Paragraph("<bullet>&bull;</bullet>"+txt, bullet_style))

    def add_spacer(self, height_in_inches):
        self.elements.append(Spacer(0.1*inch, height_in_inches*inch))

    def switch_page_template(self, template_id):
        self.elements.append(NextPageTemplate(template_id))

    def add_chapter(self, text=None):
        if not text:
            text = 'Chapter %d' % self.chapter_num
        self.add_paragraph(text, heading1_style)
        self.add_spacer(0.1)
        self.chapter_num += 1

    def build(self):
        self.doc.multiBuild(
            self.elements, canvasmaker=Canvas)



class PeptagramDoc(ReportLabDoc):
    def __init__(
            self, pdf, title, author, n_peak, mz_delta, metadata_keys):
        ReportLabDoc.__init__(self, pdf)
        self.n_peak = n_peak
        self.mz_delta = mz_delta
        self.metadata_keys = metadata_keys

        self.labels = [
            'b(3+)', 'b(2+)', 'b', 
            '>', 'aa', '<', 
            'y', 'y(2+)', 'y(3+)']

        # make title page
        self.add_spacer(4)
        self.add_paragraph(title, title_style)
        self.add_spacer(0.4)
        self.add_paragraph(author, subtitle_style)
        self.add_spacer(0.5)
        self.add_bullet("""\
            The visualizations were generated with 
            <font color="blue"><link href="http://boscoh.github.io/peptagram">peptagram</link></font>
            on %s""" % datetime.date.today())
        self.add_bullet("""\
            This is a legacy format designed to satisfy journal requirements. 
            There are better ways provided by 
            <font color="blue"><link href="http://boscoh.github.io/peptagram">peptagram</link></font>
            to visualize spectra.
            """)
        self.add_bullet("""\
            In the Peptide-Spectrum-Matches, the matched peaks were tested against the
            top %d peaks using delta(M/Z)=%f.
            """ % (self.n_peak, self.mz_delta))
        self.add_bullet("""\
            The matched peaks are identified internally in
            <font color="blue"><link href="http://boscoh.github.io/peptagram">peptagram</link></font>, 
            and differs occasionaly from the peaks matched in the search-engine. 
            """)
        self.switch_page_template('OneCol')
        self.add_page_break()

    def ion_color(self, ion_type):
        if ion_type.startswith('y'):
            return colors.blue
        elif ion_type.startswith('b'):
            return colors.red
        else:
            return colors.grey

    def add_graph(self, match):
        spectrum = match['spectrum']
        sequence = match['sequence']

        x_max = max(map(lambda x: x[0], spectrum))
        x_lims = [0, x_max*1.2]
        y_max = max(map(lambda x: x[1], spectrum))
        y_lims = [0, y_max*1.2]

        full_width = 450
        full_height = 120

        box_offset_i = 30
        box_offset_j = 30
        box_width = full_width - box_offset_i
        box_height = full_height - box_offset_j

        x_width = x_lims[1] - x_lims[0]
        y_height = y_lims[1] - y_lims[0]

        border_color = colors.black
        background_color = colors.Color(0.98, 0.98, 0.98)
        unmatched_color = colors.Color(0.8, 0.8, 0.8)

        def x_to_i(x):
            fraction = (x - x_lims[0])/float(x_width)
            return fraction*box_width + box_offset_i

        def y_to_j(y):
            fraction = (y - y_lims[0])/float(y_height)
            return fraction*box_height + box_offset_j

        drawing = Drawing(full_width, full_height)

        drawing.add(
            Rect(
                box_offset_i, box_offset_j, 
                box_width, box_height, 
                fillColor=background_color,
                strokeColor=None))
      
        for peak in spectrum[:self.n_peak]:
            i = x_to_i(peak[0])
            j = y_to_j(peak[1])
            drawing.add(
                Line(
                    i, box_offset_j, i, j, strokeColor=unmatched_color, 
                    strokeWidth=0.75))

        modifications = match['modifications']
        for ion_type in self.labels:
            if ion_type in ['<', '>', 'aa']:
                continue
            color = self.ion_color(ion_type)
            for matched_peak in self.matched_peaks[ion_type]:
                i = x_to_i(matched_peak[0])
                j = y_to_j(matched_peak[1])
                drawing.add(
                    Line(
                        i, box_offset_j, i, j, 
                        strokeColor=color, strokeWidth=0.75))
                group = Group(
                    String(
                        0, 0, matched_peak[2], fontSize=5, fontName='Helvetica',     
                        textAnchor='start', fillColor=color))
                group.translate(i+1, j+2)
                group.rotate(90)
                drawing.add(group)

        # draw labels
        drawing.add(
            String(
                box_offset_i + 0.5*box_width, box_offset_j - 10, 
                'M/Z', 
                fontSize=6, fontName='Helvetica',     
                textAnchor='middle', fillColor=unmatched_color))

        drawing.add(
            String(
                box_offset_i, box_offset_j - 10, 
                str(x_lims[0]), 
                fontSize=6, fontName='Helvetica',     
                textAnchor='middle', fillColor=unmatched_color))

        drawing.add(
            String(
                box_offset_i + box_width, box_offset_j - 10, 
                str(int(x_lims[1])), 
                fontSize=6, fontName='Helvetica',     
                textAnchor='middle', fillColor=unmatched_color))

        drawing.add(
            String(
                box_offset_i - 5, box_offset_j - 2, 
                str(y_lims[0]), 
                fontSize=6, fontName='Helvetica',     
                textAnchor='end', fillColor=unmatched_color))

        drawing.add(
            String(
                box_offset_i - 5, box_offset_j + box_height - 2, 
                str(int(y_lims[1])), 
                fontSize=6, fontName='Helvetica',     
                textAnchor='end', fillColor=unmatched_color))

        self.elements.append(drawing)

    def add_table(self, match):
        sequence = match['sequence']
        n_residue = len(sequence)

        masses = {}
        spectrum = match['spectrum']
        modifications = match['modifications']

        for ion_type in self.labels:
            if ion_type == 'aa':
                masses[ion_type] = ['aa'] + list(sequence)
            elif ion_type == '>':
                masses[ion_type] = ['>'] + map(str, range(1, n_residue)) + ['']
            elif ion_type == '<':
                masses[ion_type] = ['<',''] + map(str, range(n_residue-1, 0, -1))
            else:
                ion_list = [ion_type] + ['' for i in range(n_residue)]
                for matched_peak in self.matched_peaks[ion_type]:
                    ion = matched_peak[2].split('(')[0]
                    m = re.search(r'\d+', ion)
                    if m:
                        resnum = int(m.group(0))
                        if ion_type.startswith('y'):
                            resnum = n_residue-resnum+1
                        ion_list[resnum] = " %.1f" % float(matched_peak[0])
                masses[ion_type] = ion_list

        data = []
        for j in range(n_residue+1):
            data.append([])
            for i_label, label in enumerate(self.labels):
                data[j].append(masses[label][j])
        table = Table(data, hAlign='CENTER')

        styles = []
        border_color = colors.Color(0.9, 0.9, 0.9)
        mid_color = colors.Color(0.4, 0.4, 0.4)
        for i_label, label in enumerate(self.labels):
            if label in ['>', '<', 'aa']:
                styles.append(('ALIGN', (i_label,1), (i_label,-1), 'CENTER'))
                styles.append(
                    ('TEXTCOLOR', (i_label,0), (i_label,-1), mid_color))
            else:
                color = self.ion_color(label)
                styles.append(
                    ('TEXTCOLOR', (i_label,0), (i_label,-1), color))
                styles.append(
                    ('ALIGN', (i_label,1), (i_label,-1), 'RIGHT'))
        styles.extend([
            ('ALIGN', (0,0),(-1,0), 'CENTER'),
            ('FONTSIZE',(0,0),(-1,-1),6),
            ('LEADING',(0,0),(-1,-1),8),
            ('TOPPADDING',(0,0),(-1,-1),0),
            ('BOTTOMPADDING',(0,0),(-1,-1),0),
            ('INNERGRID', (0,0), (-1,-1), 0.25, border_color),
            ('BOX', (0,0), (-1,-1), 0.25, border_color),
        ])
        table.setStyle(TableStyle(styles))

        self.elements.append(table)
        
    def add_match_page(self, match):
        sequence = match['sequence']
        spectrum = match['spectrum']
        modifications = match['modifications']
        sequence = match['sequence']
        scan_id = match['attr']['scan_id']

        title = 'Scan %s: %s' % (scan_id, sequence)
        print('Building - ' + title)
        self.add_chapter(title)

        self.matched_peaks = {}
        for ion_type in self.labels:
            if ion_type not in ['aa', '>', '<']:
                self.matched_peaks[ion_type] = \
                    peptidemass.map_matched_ions(
                        ion_type, 
                        sequence, 
                        spectrum[:self.n_peak], 
                        mz_delta=self.mz_delta,
                        modified_aa_masses=modifications)

        # write out metadata
        if 'modifications' in match:
            modifications = match['modifications']
            mod_strs = []
            for m in modifications:
              i = m['i']
              s = "mass(%s%d)=%.1f  " % (sequence[i], i+1, float(m['mass']))
              mod_strs.append(s)
            match['attr']['modifications'] = ', '.join(mod_strs)
        for key in self.metadata_keys:  
            if key in match['attr']:
                val = match['attr'][key]
                if key == 'q_value':
                    val = "%.4f" % float(val)
                self.add_paragraph('%s: %s' % (key, val))

        self.add_spacer(0.2)

        self.add_graph(match)

        self.add_table(match)

        self.add_page_break()




if __name__=="__main__":
    match = \
        {'modifications': [{'i': 7, 'mass': 160.0306}], 'attr': {'m/z': '627.7907', 'scan_id': 2794, 'mass': '1253.5668', 'source': '2010-12-05_Yeast_Trypsin_FT_HCD_Rep3.mzML', 'q_value': 0.0358382813104469, 'morpheus_score': '9.1669', 'retention_time': '23.2851', 'modified_sequence': 'DSAVNAVC[carbamidomethylation of C]YGAK', 'mass_diff': '-0.0041'}, 'sequence': 'DSAVNAVCYGAK', 'i': 272, 'mask': 1, 'spectrum': [(185.05540466308594, 42975.28125), (256.09271240234375, 25012.060546875), (147.11276245117188, 19624.671875), (129.10215759277344, 16190.02734375), (285.15533447265625, 15223.490234375), (136.07568359375, 13320.6318359375), (381.1232604980469, 13042.6796875), (175.07125854492188, 12826.455078125), (275.17083740234375, 12584.76953125), (598.265869140625, 12018.7373046875), (203.0659637451172, 11364.2998046875), (130.08612060546875, 10944.3671875), (438.2334899902344, 10772.6748046875), (183.1492462158203, 10437.4208984375), (296.1059265136719, 9972.052734375), (186.08714294433594, 8930.7890625), (357.2502746582031, 7618.14208984375), (484.2616271972656, 7540.4755859375), (143.11769104003906, 6948.3505859375), (133.04281616210938, 6891.43603515625), (158.09201049804688, 6716.30419921875), (966.5445556640625, 6547.98486328125), (274.1029968261719, 6121.45068359375), (218.1499786376953, 5493.49365234375), (171.11276245117188, 5130.3662109375), (211.14340209960938, 5094.90185546875), (695.3556518554688, 5094.72509765625), (260.19622802734375, 4948.6103515625), (214.11831665039062, 4834.89794921875), (324.10205078125, 4703.16162109375), (110.0711441040039, 4681.60595703125), (212.10305786132812, 4477.79345703125), (260.1063232421875, 4249.50830078125), (240.13357543945312, 4194.55908203125), (130.0498046875, 4052.209228515625), (169.13357543945312, 3980.943115234375), (228.0985565185547, 3511.806884765625), (655.377197265625, 3363.16064453125), (967.555419921875, 3143.34375), (697.3330078125, 3023.61865234375), (187.1069793701172, 2923.6943359375), (583.2449340820312, 2901.6982421875), (439.23577880859375, 2863.4404296875), (881.4544677734375, 2809.265380859375), (599.2688598632812, 2796.46630859375), (226.1172637939453, 2733.450439453125), (261.156005859375, 2655.503662109375), (157.0609893798828, 2559.634765625), (244.16632080078125, 2436.15576171875), (452.16180419921875, 2435.171630859375)], 'intensity': 0.942658749903285}

    doc = PeptagramDoc(
        'out.pdf',
        "Morpheus Search Results",
        "Craig Venger",
        "test_set",
        50,
        0.1,
        ['modified_sequence', 'modifications', 'source', 
            'scan_id', 'retention_time', 'm/z', 'mass', 
            'mass_diff', 'morpheus_score', 'q_value'])
    doc.add_match_page(match)
    doc.build()
    os.system('open out.pdf')




