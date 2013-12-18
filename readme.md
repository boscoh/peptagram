# PEPTAGRAM 

`peptagram` generates a single-page HTML5 web-app to visualize proteomics analyses for:

1. an single proteomics experiment ([example][example1]),
2. a comparison between multiple proteomics experiment ([example][example2]).

[example1]:http://boscoh.github.io/peptagram/examples/single/index.html
[example2]:http://boscoh.github.io/peptagram/examples/multiple/index.html

The resulting HTML5 web-app is:

- instantly available: *runs on any modern webbrowser*
- cross-platform: *everyone has a webbrowser*
- self-contained: *just zip directory and send by email*
- publishable: *installable on any server, even on Dropbox*

User guide to the visualisation webapp at <http://boscoh.github.io/peptagram/vizhelp.html>.

`peptagram` is not meant to replace heavy-duty tools such as peptide-shaker or PepXMLViewer.  Instead, it provides a useful overview of proteomics experiments in a convenient format for end-users.

*This info in friendlier form at <http://boscoh.github.io/peptagram>.*

## Installation

`peptagram` consists of a set of python scripts that converts proteomics data into an HTML5 visualisation. 

To generate the visualisations, you must have python installed. If you have the python installer pip, then:

    > pip install peptagram

Otherwise, download and unzip the package from <https://github.com/boscoh/peptagram/archive/master.zip>

And install from the package directory:

    > python setup.py install

`peptagram` has two python dependencies: 

  1. [pymzml](https://github.com/pymzml/pymzML) to read `.mzML` files
  2. [uniprot](https://github.com/boscoh/uniprot) to get protein sequences from <http://uniprot.org>. 

These should be installed automatically from the above scripts, but if that fails, you may need to install them manually.

## Generating visualisations

`peptagram` consists of Python modules that convert proteomics data into an HTML5 visualisation. 

A set of examples can be found in this ~100 MB download: <http://monash.edu/proteomics/peptagram/examples-0.1a.zip>.

For a fully realized visualisation, `peptagram` requires:

1. peptide-spectrum matches
2. protein groupings
3. sequences of identified proteins
4. tandem MS/MS spectra 

Given the variety of proteomics formats, different combinations of data are needed to generate visualisation.

In this download, we have examples of:

- Morpheus & mzML
- TPP & fasta/uniprot
- TPP & x!Tandem
- TPP & Mascot & fasta
- MaxQuant & fasta

These will be discussed in detail below.


### Example: Morpheus & mzML

Morpheus is a search engine designed for high-quality data, where the assumption that MS/MS peaks are well-resolved results in better performance. As Morpheus does not come with a bundled viewer, `peptagram` provides a unique tool to view Morpheus proteomics search results. 

In the unzipped examples directory, the morpheus script is `run_morpheus.py`. The script first imports the relevant `peptagram` modules:

    import peptagram.morpheus
    import peptagram.mzml
    import peptagram.proteins

The following function reads in the sequences and protein groupings from `.protein_group.tsv`, the peptide-spectrum matches from `.PSMs.tsv` and the modifications from `modifications.tsv`, and put into a python data-structure called `proteins`:

    proteins = peptagram.morpheus.get_proteins(
        'morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
        'morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
        'morpheus/modifications.tsv'
        )

__OPTIONAL:__ If you did the Morpheus search with .mzML files, you can load in the raw spectra into `proteins` from the `.mzML` file, which will be displayed. The function will match the spectra in the `.mzML` to  entries in the `.PSM.tsv` file. For this case, only one data-set has been read into `proteins`, so the second parameter is given 0, to refer to the first data-set:

    i_source = 0
    peptagram.mzml.load_mzml(
        proteins, i_source, 
        'morpheus/OK20130822_MPProtomap_KO1.mzML')

Now that `proteins` is properly read, we can generate the visualation. The visualisation requires a little bit information to provide useful annotations for the user.  This is to be put in a `data` dictionary:

    data = {
      'title': 'Morpheus Example', 
      'proteins': proteins,
      'source_labels': [''],
      'color_names': ['', '', ''],
      'mask_labels': [],
    }

A description of the fields:
  
  - `title`: the text that will be displayed at the top of the web-app
  - `proteins`: the proteins data structure referred to above
  - `source_labels`: in comparison mode, this is the label for the different proteomics experiment. ignored if it is an empty list
  - `color_names`: in comparison mode, the labels for the colors in the legend
  - `mask_labels`: alternative masking for the display of higher accuracies

The `data` dictionary is then passed into the function that generates the visualisation web-app in the indicated directory:

    peptagram.proteins.make_overview_visualisation(
        data, 'morpheus/overview')

Which is run on the command-line in the examples directory:

    python run_morpheus.py


### Example: TPP & fasta/uniprot

The Transatlantic Protein Pipeline (TPP) represents one of the largest open-source proteomics toolkits. As the TPP have pushed for their `.protXML` and `.pepXML` formats as standards, the search results can come from any number of search-engines. As such, by handling TPP files, `peptagram` can  handle search data indirectly from many different sources.

The example script is called `run_tpp.py`, and the requiredn modules are:

    import peptagram.tpp
    import peptagram.proteins

First, the protein-groups and peptide-spectrum matches are read in from a `.protXML` file and potentially several `.pepXML` files, which is treated as a list of filenames:

    protein_error = 0.01

    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'tpp/hca-lysate-16.prot.xml', 
        ['tpp/hca-lysate-16.pep.xml'],
        peptide_error=0.05,
        protein_error=protein_error)

Two extra parameters are provided, `peptide_error` is a false-positive-error cutoff for the peptide entries in the `.pepXML` file, and `protein_error` is a false-positive-error cutoff for the `.protXML` file. These are calculated from the probability-error distributions included in every `.pepXML` and `.protXML` file.

As TPP files lack sequences, we must read them in externally. From our experience, sequence identifiers for protein sequences are oftern formatted differently by different programs, even if it is in essence, the same identifier. As such, there are  hooks in `peptagram` to clean up sequence-identifier for matching purposes. This is basically a idempotent string transformation function, such as:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

Now we load the protein sequences from a `.fasta` file using `clean_seqid` to transform sequence identifiers and match the identifiers from `proteins` and the `.fasta` file:

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 
        'tpp/HUMAN.fasta', 
        clean_seqid=clean_seqid,
        iso_leu_isomerism=False)

Notice the `iso_leu_isomerism` flag, which is required for some packages that don't distinguish this mass isomerism in peptides.

**Alternatively**: you can load the sequences directly from <http://uniprot.org>. This can be seen in the `run_tpp_uniprot.py` script. The sequence loading occurs with the alternative function, which requires a `cache_basename` - a filename location for temporary files:

    peptagram.proteins.load_sequences_from_uniprot(
        proteins, 
        clean_seqid=clean_seqid,
        cache_basename='tpp/uniprot')
        
The sequences are loaded into the data structure, and the positions relative to the full protein of the peptides in the peptide-spectrum-matches are calculated.

And finally, the step to generate the visualisations:

    data = {
      'title': 'TPP example',
      'proteins': proteins,
      'source_labels':source_names,
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': [protein_error],
    }
    peptagram.proteins.make_overview_visualisation(
      data, 'tpp/overview')


### Example: X!Tandem & TPP

The default search-engine that comes with the TPP is X!Tandem. The great thing about `.tandem` files is that they provide both protein sequences and the MS/MS spectra for each peptide-spectrum match. You can thus generate a full `peptagram` visualization with `.tandem`, `.protXML` and `.pepXML` files. 

The example script is `run_xtandem.py`. First, read in the `proteins` data structure from the `.pepXML` and `.protXML` files:

    import peptagram.parse
    import peptagram.xtandem
    import peptagram.tpp
    import peptagram.proteins

    errors = [0.05, 0.025, 0.01]
    
    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'xtandem/interact.prot.xml', 
        ['xtandem/interact.pep.xml'], 
        peptide_error=max(errors))

  Here, `.pepXML` was generated from 3 different `.tandem` files:

      tandems = [
        'xtandem/Seq23282_E1O1.tandem',
        'xtandem/Seq23283_E1O1.tandem',
        'xtandem/Seq23284_E1O1.tandem',
      ]

 As such, we need to make use of the `source_names` return variable from above, which contains the names of the source files in the `.pepXML` entries. Hopefully, these will match the `.tandem` to the `source_names`. To figure out the match of the `.tandem` files to the `source_names`:

    def get_i_source(tandem, source_names):
      basename = os.path.splitext(os.path.basename(tandem))[0]
      for i, source_name in enumerate(source_names):
        if basename in source_name:
          return i
      raise IOError('Couldn\'t match {} to {}'.format(basename, source_names))

Once matched, we load the `.tandem` files into `proteins`:

    for tandem in tandems:
      i_source = get_i_source(tandem, source_names)
      peptagram.xtandem.load_xtandem_into_proteins(proteins, tandem, i_source)

Then we can generate the web-app:

    data = {
      'title': 'X!Tandem example',
      'proteins': proteins,
      'source_labels': map(peptagram.parse.basename, source_names),
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': errors,
    }
    peptagram.proteins.make_comparison_visualisation(
        data, 'xtandem/comparison')


### Example: TPP & Mascot & fasta

Despite the weirdness of `.dat` files, Mascot is one of the oldest and most popular search engines. Here, we have a parser for Mascot `.dat` files, which is a mish-mash of mime-type, xml, and random acts of text. As Mascot does not group proteins, we only provide an example where mascot `.dat` files have been processed by the TPP. 

In this instance, the `.pepXML` file was generated from 2 mascot `.dat` files:

    mascot_dats = [
      'mascot/F022043.dat',
      'mascot/F022045.dat',
    ]

 So first we read in the `proteins` data structure from the `.pepXML` and `.protXML` files:

    import peptagram.parse
    import peptagram.mascot
    import peptagram.tpp
    import peptagram.proteins

    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'mascot/interact.prot.xml', 
        ['mascot/interact.pep.xml'], 
        peptide_error=None,
        protein_error=None)

We define a function to match the `.dat` files to the `source_names` in the original read with this function:

    def get_i_source(fname, source_names):
      basename = os.path.splitext(os.path.basename(fname))[0]
      for i, source_name in enumerate(source_names):
        if basename in source_name:
          return i
      raise IOError('Couldn\'t match {} to {}'.format(basename, source_names))

Then we load the MS/MS spectra from the `.dat` files:

    for mascot_dat in mascot_dats:
      i_source = get_i_source(mascot_dat, source_names)
      peptagram.mascot.load_mascot_dat_to_proteins(proteins, i_source, mascot_dat)

At this point `proteins` contain protein groupings, peptide-spectrum matches and raw spectra. We still need the protein sequences. Then we load the sequences in from a `.fasta` database with a sequence-identifier normalizng fucntion `clean_seqid`:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 'mascot/HUMAN.fasta', clean_seqid)

Finally, we generate the webapp:

    data = {
      'title': 'Mascot example',
      'proteins': proteins,
      'source_labels': map(peptagram.parse.basename, source_names),
      'color_names': ['P=1', 'P=0', ''],
      'mask_labels': [],
    }
    peptagram.proteins.make_comparison_visualisation(
        data, 'mascot/comparison')


### Example: MaxQuant & fasta

We have used MaxQuant mainly for its ability to do isotype-labelling comparisons (eg SILAC). Maxquant results are stored as tab-separated-value files as `.txt` files in a summary directory. 

The script is `run_maxquant` that reads the relevant `.txt` files from the summary directory:

    import peptagram.maxquant
    import peptagram.proteins

    proteins, sources = peptagram.maxquant.get_proteins_and_sources(
        'maxquant/summary')

A word of warning: the spectra read in is only a list of matched peaks, which is not really that useful for evaluating the quality of the peak matching. Unfortunately, the MaxQuant output was quite confusing and it was not clear how to calculate the theoretical heavy and light isotope peaks.

Then, we use the next function to use the isotype intensity ratios to populate the `intensity` field in the peptide-spectrum matches, which is used to generate the color:

    peptagram.maxquant.calculate_ratio_intensities(proteins, max_ratio=1.5)

As MaxQuant files do not contain sequences, here's the `.fasta` loading:

    def clean_seqid(seqid):
      if '|' in seqid:
        return seqid.split('|')[1]
      else:
        return seqid

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, 
        'maxquant/yeast_orf_trans_all_05-Jan-2010.fasta', 
        clean_seqid=clean_seqid)

And now we can generate the webapp, where the colors represents the SILAC ratios. As such we give the ratios as the `color_names`:

    data = {  
      'title': 'Maxquant example', 
      'proteins': proteins,
      'source_labels': sources,
      'color_names': ['1.5', '1', '0.66'],
      'mask_labels': [],
    }
    peptagram.proteins.make_comparison_visualisation(
      data, 'maxquant/comparison')


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

Loading protein sequences in the `proteins` is a necessary step in `peptogram` as it is required to generate the visualisations. As well, it is necessary to view the protein. Several of the proteomics data formats do not provide this information and so we need to read in the fasta sequences from another source.

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

