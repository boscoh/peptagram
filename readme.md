# PEPTAGRAM

`peptagram` provides a hassle-free visual overview of proteomics data for casual/non-expert users. From your proteomics search results, `peptagram` generates a single-page HTML5 web-app to visualize:

1. an entire proteomics experiment ([example]()),
2. a comparison between multiple proteomics experiment ([example]()).

The resulting HTML5 web-app is:

- instantly available: *runs on any modern webbrowser*
- cross-platform: *everyone has a webbrowser*
- self-contained: *just zip directory and send by email*
- publishable: *installable on any server, even on Dropbox*

`peptagram` is not meant to replace heavy-duty proteomics visulaization tools such as peptide-shaker or PepXMLViewer. Instead, it provides a useful overview of proteomics experiments in a package that is convenient for down-stream purposes.

As well, `peptagram` provides a number of parsers for proteomics search results.

## Installation

`peptagram` consists of a set of Python scripts that parse proteomics data and sets up the web-app. The web-app itself is a custom Javascript/HTML5 app that is included in the package. 

To generate the visualizations, you must work with Python on the commandline. You need to have installed Python and pip, the Python package installer. Then:

    pip install peptagram

`peptagram` has one dependecy, the wonderful [`pymzml`](https://github.com/pymzml/pymzML) package to read `.mzML` files.

## Generating visualizations

In proteomics, as you may very well know, there are many different formats, each choosing to describe sets of information. `peptagram` currently provide parsers to work with:

- mzML
- Morpheus
- TPP
- X!Tandem
- Mascot
- MaxQuant
- fasta 

For a fully realized visualization, `peptagram` requires:

1. peptide-spectrum matches
2. protein groupings
3. sequences of identified proteins
4. tandem MS/MS spectra 

As different formats contain different subsets of information, for a functional `peptagram` visualiation, different data files will be needed. In the following, we will describe a number of cases we have got working.

### EXAMPLE: Morpheus and mzML

`morpheus` is a search engine designed for high-quality data, where the assumption that MS/MS peaks are well-resolved results in much better performance. As `morpheus` does not come with a bundled viewer, `peptagram` provides a unique tool to view `morpheus` results. 
First create a Python script in a text editor with a name like `morpheus_example.py`. This will read a set of `morpheus` protein groupings and match the proteins with peptide-spectrum matches. Make sure the filenames are correct relative to the directory where you save `morpheus_example.py`. 

`morpheus` generates protein groupings and the nice thing about `morpheus` is that protein sequences and descriptions are included in the `.protein_group.tsv` files. This means the primary data structure `proteins` can be generated with the `morpheus` `.protein_groups.tsv` and `.PSMs.tsv` files. Also required is the `modifications.tsv` file:

    import peptagram.morpheus
    import peptagram.mzml
    import peptagram.proteins

    proteins = peptagram.morpheus.get_proteins(
        'example/morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
        'example/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
        'example/morpheus/modifications.tsv'
        )


__OPTIONAL:__ If you did the `morpheus` calculation with `.mzML` files, then there is an optional step where you can load in raw spectra from the `.mzML` file. The visualization will still work without this step, but no spectra will be displayed:

    peptagram.mzml.load_mzml(
        proteins, 0, 'example/morpheus/OK20130822_MPProtomap_KO1.mzML')

Now that we have loaded the peptide-spectrum matches and protein groups into  `proteins`, we can generate the web-app visualation. This requires a population of a `data` dictionary that contains parameters that will be passed into the web-app:

    out_dir = 'out/morpheus-pr'
    data = {
      'title': 'Morpheus Example', 
      'proteins': proteins,
      'source_labels': [''],
      'color_names': ['1.0', 'score/n', ''],
      'mask_labels': [],
    }
    peptagram.proteins.make_proteins_directory(data, out_dir)

A quick description of the fields:
  
  - title: the text that will be displayed at the top of the web-app
  - proteins: the proteins data structure, in the future you may want to construct your own
  - source_labels: in comparison mode, this is the label for the different proteomics experiment. ignored if it is an empty list
  - color_names: in comparison mode, the labels for the colors in the legend
  - mask_labels: alternative masking for the display of higher accuracies

On the command-line, run the script to generate the visualization:

    python morpheus_example.py

### EXAMPLE: TPP with fasta and mzML 

The Transatlantic Protein Pipeline (TPP) represents one of the largest open-source proteomics toolkits. As the TPP have pushed for their `.protXML` and `.pepXML` formats as standards. `.protXML` and `.pepXML` search results can come from any number of search-engines. `peptagram` can generate visualizations from these files if the protein sequences are also available in the form of `.fasta` files:

    import peptagram.tpp
    import peptagram.mzml
    import peptagram.fasta
    import peptagram.proteins

    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'tpp/hca-lysate-16.prot.xml', 
        ['tpp/hca-lysate-16.pep.xml'],
        peptide_error=0.05,
        protein_error=0.01)

However, for our purposes, these files are lacking the protein sequences required for the visualization, so protein sequences must be loaded in from `.fasta` files. Before we do that, we will define a sequence ID function:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

The reason for this is simply that noone seems to be agree on the exact format of sequence identifiers, and so the pepXML and fasta files may differ ever so slightly in format, thereby rendering the matching impossible. As such, a seqid cleaning function is always available to make the seqids consistent, in both directions. Then 

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 
        'tpp/HUMAN.fasta', 
        clean_seqid=clean_seqid,
        iso_leu_isomerism=False)

__OPTIONAL:__ As before, if `.mzML` files were used to generate the search results, then the MS/MS spectra can be read in to generate the spectra visulaizations:

    peptagram.mzml.load_mzml(
        proteins, 0, 'example/morpheus/OK20130822_MPProtomap_KO1.mzML')

And finally, the step to generate the visualiations:

    data = {
      'title': 'TPP example',
      'proteins': proteins,
      'source_labels': ['hca'],
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': map(str, errors),
    }
    peptagram.proteins.make_proteins_directory(data, 'out/tpp-pr')

On the command-line, run the script to generate the visualization:

    python tpp_mzml_example.py


### EXAMPLE: X!Tandem in TPP

The default search-engine that comes with the TPP is X!Tandem. It is easy to generate  visualizations with `.protXML` and `.pepXML` files that have been generated from `.tandem` files. This is because `.tandem` files contain protein sequences and MS/MS peaks.

    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'xtandem/interact.prot.xml', 
        ['xtandem/interact.pep.xml'], 
        peptide_error=max(errors))

In this particular example, the pepXML was generated from 3 different .tandem files. To read this into the `proteins` data structure, we will need to match the .tandem file to the `source_names` returned by the `get_proteins_and_sources` function. So here's a function to figure out the match between the `tandem` file and `source_names`

    def get_i_source(tandem, source_names):
      basename = os.path.splitext(os.path.basename(tandem))[0]
      for i, source_name in enumerate(source_names):
        if basename in source_name:
          return i
      raise IOError('Couldn\'t match {} to {}'.format(basename, source_names))

So, here are the three original tandems:

      tandems = [
        'xtandem/Seq23282_E1O1.tandem',
        'xtandem/Seq23283_E1O1.tandem',
        'xtandem/Seq23284_E1O1.tandem',
      ]

Once matched, we can load .tandem into the `proteins` data structure:

    for tandem in tandems:
      i_source = get_i_source(tandem, source_names)
      peptagram.xtandem.load_xtandem_into_proteins(proteins, tandem, i_source)

### EXAMPLE: TPP with fasta from Mascot

Mascot is one of the oldest search engines and has gone through considerable changes over the years. Here, we have a parser for the Mascot `.dat` format, which is miss-mash of mime-type, xml, and random acts of text. Mascot however, does not group proteins, and mascot `.dat` files are often run through the TPP. Here, `peptagram` can read TPP-generated `.pepXML` and `.protXMl` files, and extract the MS/MS peaks from the mascot `.dat` files. Furthermore, the protein sequences must be read from a .fasta file.

So first we read in the `proteins` data structure from the .pepXML and .protXML files:

    import peptagram.parse
    import peptagram.mascot
    import peptagram.tpp
    import peptagram.proteins


    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'mascot/interact.prot.xml', 
        ['mascot/interact.pep.xml'], 
        peptide_error=None,
        protein_error=None)

In this instance, the .pepXML file is generated from 2 mascot .dat files:

    mascot_dats = [
      'mascot/F022043.dat',
      'mascot/F022045.dat',
    ]

We define a function to match the `mascot_dat` files to the `source_names` in the original read with this function:

    def get_i_source(fname, source_names):
      basename = os.path.splitext(os.path.basename(fname))[0]
      for i, source_name in enumerate(source_names):
        if basename in source_name:
          return i
      raise IOError('Couldn\'t match {} to {}'.format(basename, source_names))

Then we load the MS/MS spectra from the .dat files:

    for mascot_dat in mascot_dats:
      i_source = get_i_source(mascot_dat, source_names)
      peptagram.mascot.load_mascot_dat_to_proteins(proteins, i_source, mascot_dat)

Finally, we define a `clean_seqid` function to handle sequence identifiers:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

then we load the sequences in:

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 'mascot/HUMAN.fasta', clean_seqid)

Generate the webapp:

    data = {
      'title': 'Mascot example',
      'proteins': proteins,
      'source_labels': map(peptagram.parse.basename, source_names),
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': ['1.0'],
    }
    peptagram.proteins.make_peptograph_directory(
        data, 
        'mascot/webapp')


### EXAMPLE: MaxQuant

We have used MaxQuant mainly for its ability to do isotype-labelling comparisons (SILAC for instance). Maxquant results are stored as tab-separated-value files as `.txt` files in a summary directory. This contains protein groups, peptide-spectrum matches, and a list of matched peaks. Nevertheless, this list of matched peaks is not useful as it does not give a visual representation of the fit of the data to the raw spectra.

Nevertheless, this can be read as: 

    import peptagram.maxquant
    import peptagram.mzml
    import peptagram.fasta
    import peptagram.proteins

    proteins, sources = peptagram.maxquant.get_proteins_and_sources(
        'maxquant/summary')

To extract the ratios, and use them for coloring:

    peptagram.maxquant.calculate_ratio_intensities(proteins, max_ratio=1.5)

However, MaxQuant files do not contain sequences, and so:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 
        'maxquant/yeast_orf_trans_all_05-Jan-2010.fasta', 
        clean_seqid=clean_seqid)


## User Guide for Visualizations

The visalization presents an integrated display of:

- peptide distributions for identified proteins
- single-sequence/comparison-distribution display
- peptide-spectrum matches
- optional top 50 MS/MS peaks for each peptide-spectrum match

_Protein List._ All identified proteins are shown on the left, withe the peptides illustrated in a the bar graph. All proteins are normalized to the same width of the screen as there is too great a difference between protein sequence lengths to show this in proportion. Clicking on a protein on this light will bring up a detailed display of the sequence, with the peptides highlighted by clickable links. This list can be sorted by the pop-up menu at the top of the section. All attributes displayed in the protein info panel can be sorted against.

_Keyboard Shortcuts._ Several keyboard shortcuts have been included, and these are labeled above the sub-headings of the relevant section. These are:

   - N and P: to move through proteins in the protein list
   - J and K: to move through the peptide-spectrum matches
   - U and D: to move between experiments in the comparison mode

_Searching_ can be done using the built-in web-browser search. All the text for the proteins, and principal seqid are always available in the protein-list view, and thus searchable.

_Peptide-Spectrum Matches._ Below that is the list of peptide spectrum matches, listed by their positions in the sequence and the sequnce. If modifications were read in, they will be highlighted. Click on thes will display the MS/MS spectra if available.

_Spectrum Viewer._ If the MS/MS spectra have been loaded then, if a peptide-spectrum match is clicked, a simple spectrum viewer is shown of the top 50 peaks found in the peptide-spectrum match. Below the viewer is an interactive ion table, showing b-ions and y-ions. By click on the labels at the top of the table, different charge states for the b-ions and y-ions are displayed, both in the the table and the spectrum viewer. Clicking on the viewer will zoom to the nearest peak, where the viewer will scale to the height of the zoomed peak, which should now be in the middle of the viewer. Clicking again on the viewer will zoom out. The viewer automatically scales to the minimum-maxium peak m/z values.

### Single Experiment Mode

In the single-experiment mode, the central panel shows the entire sequence of the protein rendered in an easy to scan 5 x 5 blocks, with the the peptides identified highlighted as clickable links

### Experiment Comparison Mode

In the experiment-comparison mode, the central model shows an interactive display where, for a selected protein, the central panel shows a display consisting of rows. Each row represents the peptides identified for each separate experiment.

## Python Programming API

### Basic Data Structure: proteins

The main focus of the Python scripts is turn all the different proteomics data format into one consistent data-structure, which is a JSON-compatible dictionary called `proteins`. In YAML format, `proteins` is laid out something like this:

    'example-seqid':
      sequence: 'ACDEFGHKLMNP'
      description: 'Some tasteless protein'
      attr:
        key1: value1
        key2: value2
        other_seqids: 
          - another-seqid
        seqid: 'example-seqid'
      sources: 
        -
          peptides:
            - 
              sequence: 'CDE'
              i: 1
              intensity: 1.0
              mask: 0.0
              attr:
                key3: value3
                key4: value4
              spectrum:
                -
                  - 501
                  - 34.3
                -
                  - 503.4
                  - 82.3

This data-structure is written verbeten into a JSON-based javascript file, which is then processed by the javascript application that generates the visual display. Optional fields are moved into the `attr` dictionary.

The primary level of organization of the dictionary `proteins` is, obviously, at the level of proteins. All peptide-spectrum matches will be sorted into each protein match. Presumably, the results will have been sorted into protein groups and read in with a usable sub-set of representative proteins for the protein identification.

The peptide information are sorted into separate lists in the `sources` field. This allows a clear demarcation for different experiments. However, for single experiments, this adds an extra layer, where the peptides must be accessed as:
 
    source = proteins['example-seqid']['sources'][0]
    peptides = source['peptides']


### Peptide-spectrum match list

Each peptide entry in the `peptides` list represents a distinct peptide-spectrum match. 

    sequence: 'CDE'
    i: 1
    intensity: 1.0
    mask: 0.0
    attr:
      modifications: 
        -
          i: 0
          mass: 344.4
      key3: value3
      key4: value4
    spectrum:
      -
        - 501
        - 34.3
      -
        - 503.4
        - 82.3

`i` gives 0-based position of the peptide, and should match that of the full sequence.

`intensity` gives a value from -1.0 to 1.0, which is used to generate colors in the experiment-comparison mode. The coloring goes from a high color associated with 1.0 to a neutral color associated with 0.0, and down to a color associated with -1.0

`mask` gives an optional value to control masking for certain peptide-spectrum matches. The webapp can accept several master `mask` values. When a master `mask` value is selected, all peptide-spectrum matches with higher `mask` values are hidden.

`spectrum` gives a list peaks that will be used in the spectrum viewer. The first number corresponds to the m/z value and the second, to the intensity. In general, only the top 50 are read in, but of course you can add more, which will probably bloat the javascript file.

`modifications` is an optional field that describes any amino acid modifications in the peptide. It is a list of dictionaries. In each dictionary, the `i` gives the position, and `mass` gives the mass of the modified amino acid. To allow for N-terminal modifications, `i` can take the value -1. For C-terminal modifiactions, `i` can take the value n where n is the length of the sequence.

### Sequence identifiers 

As discussed in the Examples, handling sequence identifiers (seqids) correctly is a recurring problem in bioinformatics. If seqids are not formatted correctly, it becomes impossible to match data from different sources. So it useful if we can format or transform all seqids into a consistent format. This idea of passing in seqid transforming functions is available in all parts of `peptagram`. This way you can organize the seqids as you read them in. 

In `peptagram.proteins`, there is a useful convenience function `change_seqids_in_proteins` that transforms the seqids found in a `proteins` datastructure, including the alternate seqids found in `other_seqids`. It takes a any string function `clean_seqid` and transforms all seqids found in `proteins` to this function, including the protein seqid keys at the top level.

    import peptogram.proteins

    def change_seqid(seqid):
      return seqid.split('|')[0]

    peptagram.proteins.change_seqids_in_proteins(proteins, clean_seqid)

### Protein sequences

Loading protein sequences in the `proteins` is a necessary step in `peptogram` as it is required to generate the visualizations. As well, it is necessary to view the protein. Several of the proteomics data formats do not provide this information and so we need to read in the fasta sequences from another source.

The most common source is a `.fasta` file, preferably one that was used for the peptide search. To load the correct sequence, the sequence identifiers between the peptide-spectrum matches and the fasta files must much. We offer where the protein sequences used in the peptide serch are The different data formats AnAdd uniprot metadata, and reorganization of data.

    import peptogram.proteins

    def change_seqid(seqid):
      return seqid.split('|')[0]

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 'mascot/HUMAN.fasta', clean_seqid)

### Multiple data sources

As the structure of each protein in `proteins` contains potentially several sources, it's important to track the name of the different experiments. Most of the parsing methods returns a `source_names` list that contains the name of the different source files. Once cleaned up, these can serve as labels for the different experiment, especially in the view of the experiment-comparison mode. Often the source name is a long directory name, where the unique part is in the basename. A quick way to clean this up is:

    import os
    source_labels = [os.path.basename(name) for name in source_names]

This `source_labels` can then be placed directly in the `data` structure used to generate the web-app.    

### Finangling MS/MS spectra

Being able to easily view spectra is one of the most useful aspects of `peptagram`. It was important that we could reconcile the scan identifiers in the different data formats so that we could load the spectra directly if the `.mzML` files available. Being able to do this has expanded the flexibility of the system. By the time each parser has generated a `proteins` structure, the peptide-spectrum matches should contain a valid `.mzML` scan identifier. This is the `scan_id` filed in the `attr` dictionary of each peptide entry.

Given this, we can run:

    peptagram.mzml.load_mzml(
        proteins, 0, 'example/morpheus/OK20130822_MPProtomap_KO1.mzML')



### Loading data into the webapp

Once the `proteins` data is appropriately filled in, we can pipe it through to the web-app generating method. For display, we need a little bit more information, and this is used to populate a `data` dictionary:

    data = {
      'title': 'Mascot example',
      'proteins': proteins,
      'source_labels': map(peptagram.parse.basename, source_names),
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': ['1.0'],
    }

`title` gives the title that will be displayed across the top of the web-app

`proteins` is the dictionary that you generated consistent with the structure given above

`source_labels` is an optional list of the names of the different experiments that will be used in the experiment-comparison mode

`mask_labels` are a list of values that will be used to mask the display of the peptides. If this list empty, then no masking will be offered.

Finally, to generate the webapp for the single-experiment mode:

And for the experiment-comparison mode:

    peptagram.proteins.make_peptograph_directory(
        data, 'mascot/webapp')

