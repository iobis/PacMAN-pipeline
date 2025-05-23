# PacMAN-pipeline
## Bioinformatics pipeline for the PacMAN project

This is the bioinformatics pipeline developed for the PacMAN (Pacific Islands Marine Bioinvasions Alert Network) project. This pipeline cleans and classifies sequences from eDNA samples. The steps in this pipeline are compiled from publicly available bioinformatic pipelines like [ANACAPA](https://github.com/limey-bean/Anacapa), [tourmaline](https://github.com/lukenoaa/tourmaline), [tagseq-qiime2-snakemake](https://github.com/shu251/tagseq-qiime2-snakemake), [pema](https://github.com/hariszaf/pema), [CASCABEL](https://github.com/AlejandroAb/CASCABEL) and [MBARI-BOG](https://github.com/MBARI-BOG/BOG-Banzai-Dada2-Pipeline). The pipeline is based on the snakemake workflow management system.

The pipeline has the following steps:

  1. **fastp** - Quality trimming and removal of sequencing adapters, chosen for the possibilty to work with NovaSeq data
  2. **Cutadapt** - Removal of primers
  3. **dada2** - ASV inference
  4. **RDP classifier** - Taxonomic assignment with a bayesian based probability method
  5. **vsearch** - Search and print out high-confidence local alignments for confirmation of RDP classifier results
  6. **BLAST** - Optional: blast search of remaining unknown sequences agains the NCBI nt database
  6. **Data formatting** - Export to DwC-A compatible tables

## How to use the pipeline
### Data preparation

The following files are required to run the pipeline:

1. Manifest
   - The manifest lists the sequence files and links them to the sample identifiers
   - Contains the columns: `sample-id`, `file-path` and `direction` (forward or reverse)
   - Example: [manifest.csv](data/config_files/manifest.csv)
2. Sample metadata
   - Contains fields that should be added to the Occurrence and DNADerivedData tables
   - :warning: **Note** Control samples can be marked by adding `occurrenceStatus` as `absent`
   --> The ASVs from these samples will be removed from all samples, before the occurrence table is made
   - Example: [sample_data.csv](data/config_files/sample_data.csv)
3. Raw sequences
   - Paired-end fastq files
4. Reference databases
   - For RDB classifier
      - A trained classifier like those available here, or build by following this tutorial
      - [Download example database here](https://pacman-databases.s3.us-east-1.amazonaws.com/COI_terrimporter.zip)
   - For vsearch
      - A fasta file in the sintax format
      - [Download example database here](https://pacman-databases.s3.us-east-1.amazonaws.com/COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_sintax.fasta)
5. Run configuration
   - Modify the [config.yaml](data/config_files/config.yaml) file for the specific run
   - `PROJECT` name: Usually a specific sample set
   - `RUN` name: the run with a specific combination of samples and/or parameters for the analysis
   - `SAMPLE_SET`: manifest file path
   - `sample-data-file`: manifest file path
   - reference database: name of the database, fasta file and taxa file.
   - Primers used in both forward and reverse configuration
   - Chosen parameters for each step. (Template file configured for CO1 data using the Leray-Geller primer set).

> :rocket: Quick start - If you want to try out the pipeline with example data, use the example config files and download the example reference databases. Also download the following data files: [forward](https://pacman-databases.s3.us-east-1.amazonaws.com/USP-24-01-172_S172_L001_R1_001.fastq.gz) and [reverse](https://pacman-databases.s3.us-east-1.amazonaws.com/USP-24-01-172_S172_L001_R2_001.fastq.gz). Use the following folder structure:
> 
> ```
> └── data
>     ├── config_files
>     │   ├── config.yaml
>     │   ├── manifest.csv
>     │   └── sample_data.csv
>     ├── raw_sequences
>     │   ├── USP-24-01-172_S172_L001_R1_001.fastq.gz
>     │   └── USP-24-01-172_S172_L001_R2_001.fastq.gz
>     └── reference_databases
>         ├── COI_ncbi_1_50000_pcr_pga_taxon_derep_clean_sintax.fasta
>         └── COI_terrimporter
>             ├── bergeyTrainingTree.xml
>             ├── genus_wordConditionalProbList.txt
>             ├── logWordPrior.txt
>             ├── rRNAClassifier.properties
>             └── wordConditionalProbIndexArr.txt
> ```

### Option 1: run with snakemake

Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). The pipeline is run with the [--use-conda](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) flag, where each rule has an isolated environment, that will be installed in the working directory using conda.

Download the repository and navigate to the root directory. Add your data files and reference databases to the `data` folder in the project root. Modify the example config file, manifest, and sample metadata file. Ensure that the reference databases are properly linked to from the config file, and that the file paths in the manifest match the location of your data files.

Perform a dry-run of the pipeline for testing with:

```
snakemake --use-conda --configfile data/config_files/config.yaml --rerun-incomplete --printshellcmds --cores 1 -np
```

Removing the `-np` flag will initiate the run.

### Option 2: run with Docker

The pipeline can also be run with Docker. In this case you will not need a copy of the pipeline code. After [installing Docker](https://www.docker.com/products/docker-desktop/), prepare a `data` folder in your workspace as explained in the data preparation section above. Then run the pipeline with the following command (change the paths to your workspace folder):

```bash
docker run --platform linux/amd64 -v $(pwd)/data:/pipeline/data -v $(pwd)/results:/pipeline/results -v $(pwd)/.snakemake:/pipeline/.snakemake pieterprovoost/pacman-pipeline
```

On Windows, you may have to use `%cd%` instead of `$(pwd)`.

### Example output

Download example output for a single sample [here](https://pacman-databases.s3.us-east-1.amazonaws.com/example_dataset_coi.zip).

## About the pipeline

The snakemake workflow will run the following jobs:

### Initiate file structure

The run will first initiate a folder structure in the `results` folder as follows

```
PROJECT
├── samples
|      ├── sample_1
|      │   ├── fastp
|      │   ├── qc
|      │   └── rawdata
|      │       ├── forward (link to sample file)
|      │       └── reverse (link to sample file)
|      ├── sample_2
|      ├── sample_3
|      |     .
|      |     .
|      |     .
|      ├── sample_n
|      └── multiqc_RUN.html
└── runs
       └── RUN
           ├── 02-cutadapt
           ├── 03-dada2
           ├── 04-taxonomy
           ├── 05-dwca
           └── 06-report
```

Samples will be linked to the file structure and their quality will be analysed with fastqc. All quality files of the raw sequence files will be summarized with multiqc, and can be found in `/PROJECT/samples/multqc_RUN.html`.

### Trimming and removing primers

The sequences are trimmed and primers are removed utilising fastp and cutadapt. fastp automatically recognizes and removes illumina adapters. fastp is used as it also detects common artifacts of NovaSeq datasets. The primers must be added to the config file in both forward and reverse (reverse complement) directions.

### dada2

ASVs are inferred with dada2, which is run in 2 steps. Initially filtering of samples is done based on user-defined parameters. The quality of sequences before and after this filtering is shown in aggregate in 2 plots (`06-reports/dada2`), and can be found separately for each samples in the `03-dada2/quality` folder.

dada2 returns the ASV-table (`03-dada2/seqtab-nochim.txt`), as well as the sequences of each asv (`03-dada2/rep-seqs.fna`). In addition the number of reads filtered at each step and remaining after sample processing are returned in the `06-report/dada2_stats.txt` file, and will be added to the report.

### Taxonomy

Taxonomy is annotated using two different methods. Annotation using the naive bayesian classifier rdp is the main method, and the results are filtered based on a user-defined probability threshold (i.e. 0.8). The output provides a confidence threshold for each taxonomic level. For rdp a trained classifier is required. Many trained classifiers can be found already online (see for example https://github.com/terrimporter/CO1Classifier), but a classifier can also be built from a custom reference database by following the steps in the tutorial of [John Quensen](https://john-quensen.com/tutorials/training-the-rdp-classifier/). 

To confirm the results with a local alignment to the reference database, also vsearch is used. The parameters of vsearch can be defined by the user. We use it to find matches with very high sequence similarity to the query (97-100%), to confirm the results of rdp as well as detect any disrepancies, or unclear results in the reference database. The top 10 hits are filtered with a built-in lca algorithm, to find the consensus taxonomy. The reference database format required for this process is the sintax format.

**Note** As two different taxonomic assignment methods are used to enable confirmation of the results, two different reference database formats are also necessary.

### BLAST and LCA (optional)

There is an option in the pipeline to further classify sequences that remained unclassified with BLASTn against the full ncbi nt database. However multiple resources are required for this to work. We recommend having a local copy of the full NCBI nt database available to run this step with the pipeline. If you have in total <10kbp of data (50 unknown sequences of 200 bp), you may also run the query in remote mode. We may include a loop to do this with more data at a later stage, but running the analysis remotely for more sequences will require a lot of time.

The user will need to also have access to NCBI-nt to TaxonID mapping files to get the scientific names of the sequences. The pipeline uses [BASTA](https://github.com/timkahlke/BASTA) to filter and classify the Blast results based on an lca analysis. If a tax database is not provided in the configfile, the pipeline will prompt BASTA to download the tax_db (gb) to the `resources` folder. This will also take a long time.

### DwC-A output

In the final steps of the pipeline LSIDs are defined for the assigned taxonomic names, and the occurrence table and dna-derived data extension table are built for submitting into OBIS.

This step also returns a table `05-dwca/Taxa_not_in_worms.csv`, containing the taxonomic names and linked asvs that were not given an lsid. This table will require manual inspection, and possibly contacting the WoRMS team.

In this step the unknown sequences are given the ID for 'Incertae sedis'. Non-marine species (most taxa with no lsid), and ASVs found in the control sample(s) (marked by occurrenceStatus: absent) are not added to the final dwca-tables. All of these can still be found in the table `05-dwca/Full_tax_table_with_lsids` as well as the `05-dwca/phyloseq_object.rds`, which can be read with the phyloseq R package for further analysis and visualization.

**Note!** With this strategy, sequences that are known but not marine, are not included in the occurrence tables, while sequences that are not known are always included (as 'Incertae sedis').

Fields that also still need to be added/modified based on the genetic data guidelines are:

  - `identificationReferences`: Website of this pipeline

### Reporting

An HTML report is made in the final steps with the statistics of the full run, to give an overview of what was done during the analysis and what the effect was on the results. Still more analysis will be added to this report.
