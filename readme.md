# PEPTAGRAM

`peptagram` provides a hassle-free overview  of proteomics data for casual/non-expert users. From your proteomics search results, `peptagram` generates a single-page HTML5 web-app to visualize:

1. an entire proteomics experiment ([example]()),
2. a comparison between multiple proteomics experiment ([example]()).

As well, `peptagram` provides a number of parsers for proteomics search results.

The resulting HTML5 web-app is:

- instantly available: *runs on any modern webbrowser*
- cross-platform: *everyone has a webbrowser*
- self-contained: *just zip directory and send by email*
- publishable: *installable on any server, even on Dropbox*

`peptagram` is not meant to replace heavy-duty proteomics visulaization tools such as peptide-shaker or PepXMLViewer. Instead, it provides a useful overview of proteomics experiments in a package that is convenient for down-stream purposes.

## Installation

`peptagram` consists of a set of Python scripts that parse proteomics data and sets up the web-app. The web-app itself is a custom Javascript/HTML5 app that is included in the package. 

To generate the visualizations, you must work with Python on the commandline. You need to have installed Python and pip, the Python package installer. Then:

    pip install peptagram

`peptagram` has one dependecy, the wonderful [`pymzml`](https://github.com/pymzml/pymzML) package to read `.mzML` files.

## Generating visualizations

In proteomics, as you may very well know, there are many different formats, each choosing to describe sets of information. `peptagram` currently provide parsers to work with:

- mzML
- Morpheus
- TransatlanticProteomicsPipeline
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

First create a Python script in a text editor with a name like `morpheus_example.py`. This will read a set of `morpheus` protein groupings and match the proteins with peptide-spectrum matches. Make sure the filenames are correct relative to the directory where you save `morpheus_example.py`:

    import peptagram.morpheus
    import peptagram.mzml
    import peptagram.proteins

    proteins = peptagram.morpheus.get_proteins(
        'example/morpheus/OK20130822_MPProtomap_KO1.protein_groups.tsv',
        'example/morpheus/OK20130822_MPProtomap_KO1.PSMs.tsv',
        'example/morpheus//modifications.tsv'
        )

The nice thing about `morpheus` is that protein sequences and descriptions are included in the `.protein_group.tsv` files. This data-structure can be used to generated the visualization. 

__OPTIONAL:__ If you did the `morpheus` calculation with `.mzML` files, then there is an optional step where you can load in raw spectra from the `.mzML` file. The visualization will still work without this step, but no spectra will be displayed:

    peptagram.mzml.load_mzml(
        proteins, 0, 'example/morpheus/OK20130822_MPProtomap_KO1.mzML')

Now that we are satisfied that the peptide-spectrum matches and protein groups are loaded into the `proteins` data structure, we can generate the web-app visualation:

    out_dir = 'out/morpheus-pr'
    data = {
      'title': 'Morpheus Example', 
      'proteins': proteins,
      'source_labels': [''],
      'color_names': ['1.0', 'score/n', ''],
      'mask_labels': [],
    }
    peptagram.proteins.make_proteins_directory(data, out_dir)

On the command-line, run the script to generate the visualization:

    python morpheus_example.py

### EXAMPLE: TPP with mzML

The Transatlantic Protein Pipeline (TPP) represents one of the largest open-source proteomics toolkits. They have pushed for their `.protXML` and `.pepXML` formats as standards. As such `.protXML` and `.pepXML` search results can come from any number of search-engines and data sources.

    import peptagram.tpp
    import peptagram.mzml
    import peptagram.fasta
    import peptagram.proteins

    proteins, source_names = peptagram.tpp.get_proteins_and_sources(
        'example/tpp/hca-lysate-16.prot.xml', 
        ['example/tpp/hca-lysate-16.pep.xml'],
        peptide_error=max(errors))

However, for our purposes, these files are lacking the protein sequences required for the visualization, so protein sequences must be loaded in from `.fasta` files.

    peptagram.proteins.load_fasta_db_into_proteins(
        proteins, '../db/HUMAN.fasta')

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



Optionally 

search results in ProteinProphet prot.xml files containing all protein groups, and peptides. The peptide-spectrum matches are contained in the PeptideProphet pep.xml files. However, neither prot.xml or pep.xml contain full sequence information, so a sequence database in FASTA format is required. Optionally, the mzML corresponding to the peptide-spectrum matches can be used to integrate the MS/MS spectrum for each match.

### EXAMPLE: TPP from X!Tandem

### EXAMPLE: Mascot with TPP

### EXAMPLE: MaxQuant

## User Guide for peptagram visualizations

Probably the first thing you want to do is check out an example: example if you are using a modern browser (tested against Chrome and Safari) [ref](). This visalization presents an integrated display of:

- peptide distributions for identified proteins
- comparison across experiments
- sequence display
- peptide-spectrum matches
- top 50 MS/MS spectrum of each match

In another mode, the visualization focuses on the sequence, with the whole sequence shown in block form with the peptides highlighted clearly. 

Keyboard shortcuts are included to easily cycle through the proteins, experiments, and peptide-spectrum matches.

All the proteins are loaded into the page. To search, simply use the in-page search of the browser, which will take you to the relevant protein entry.

Sorting is carried out through the pop-up menu for proteins over protein attributes.


## Peptagram Python API


Sequence identifier nightmare.

In order to work with different formats, Python scripts convert the different data formats into a single Python data structure, which is referred to as `proteins`. The structure of `proteins` in YAML format is:

    'protein-seqid':
      sequence: 'ACDEFGHKLMNP'
      attr:
        key1: value1
        key2: value2
        other_seqids: []
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
              peaks:
                -
                  mass: 500
                  intensity: 34.3

Essentially, the Python scripts organises the peptide-spectrum matches at the level of proteins. This facilitates the analysis of the peptide-spectrum matches from a protein-identity point of view. 

It is a container of protein, with the only key fields is 'sequence' and 'attr' contains a list of variable 'attr'. There is an optional 'other_seqids' fields in 'attr' which is quite useful when matching protein groups from different reads.

Because one of the goals is to integrate peptide-spectrum matches across experiments, the peptide-spectrum matches are separated in terms of a list of `sources`. Each source is a dictionary itself, to allow future addition of information for each `source`. The key field is `peptides` that contains a list of a `peptide` dictionary.

The `peptide` dictionary contains all the information of

masks

intensities: colors

multiple sources: source_labels, sources

JSON data structure

attr dictionary can take any JSON values. For protein attr values, these can be sorted through within the page.

python script to turn search results into a python data structure.

This data structure is piped into a json structure, that is stored within a javascript structure in the website. 

Examples are given to combine the results file in order to generate the `proteins` data structure. This is then written as an equivalent JSON object in a JAVASCRIPT object  called `data` that is loaded into a template web-page.

Add uniprot metadata, and reorganization of data.

There are quite a few packages that read proteomics (pyteomics, ) search results. Why another one? Well, one thing we are focusing on here is that the results are organized in terms of consistent data structure, abstracted over all the different search engines. The focus here is on proteins, and protein groups of related proteins. We also want to express the data structure in a convenient data structure to explore within Python.
pp