# DARPAmetapipe

This code is adapted from https://github.com/mvvucak/LASV_Sample_Classification
The collection of scripts and snakemake pipelines are used to process Illumina reads from the multimammate rodent (Mastomys species) into taxonomically classified contigs. 
In this case, the pipeline was used to identify all herpesvirus, mammarenavirus, arenaviridae, herpesvirales and adenoviridae sequences present in samples.

## Dependencies

### Conda Environment

metagenomics.yml contains most of the packages needed to run the pipelines. You will need to clone the environment using Conda.

Conda can be installed from: https://docs.conda.io/en/latest/miniconda.html

Once conda has been installed, you can create the environment with:

`conda create --name metagenomics --file metagenomics.yml`

Then activate it:

`conda activate metagenomics`


### DIAMOND Database

 The contig classification step relies on a functional Diamond2 database constructed from all RefSeq Protein entries in NCBI. Additionally, the database must also include taxonomic information for each entry (TaxID). 
 Once you have the database, simply change the following pointer in the config.yml file:

 `diamond_db:
  /home2/mvv1e/Databases/RefSeq_Protein.dmnd`


## Usage

The pipeline consists of 4 Snakemake pipelines to be run sequentially. A config.yml supplies supporting information. 
The Snakemake files should be run from the project directory using:

`snakemake --snakefile {filename} -j {allocated threads}.`


### Config File

config.yml stores the following:

#### Sample/Animal Information:

- Animal IDs under `animal` (five-digit numbers with leading 0s e.g. 00019)
- Sample IDs under `all_samples`, starting with animal ID followed by sample type and Illumina run identifiers:
    - 00019_BL_L1_NA_S13 denotes blood sample from animal 00019
- Sample IDs are grouped under their corresponding animal under `all_samples_by_animal`

#### Target taxon information:

- `tax_ranks` lists the taxonomic ranks to be identified for each Diamond hit. Used by the get_tax_ranks.py script
- `target_taxons` lists information for any taxa being searched for in the de novo output: 
    - Scientific name (e.g. mammarenavirus)
    - Taxon rank (e.g. genus)
    - TaxID (e.g. 1653394)

#### Mundane processing information:

- Expected file extensions (e.g. `fastq`, `fq.gz`)
- Database paths.

### Snakemake Pipeline

 The pipeline consists of 4 separate Snakemake files, to be run sequentially. 

#### Snakereads

This Snakefile runs basic read processing steps on raw Illumina reads, including:
- Adapter trimming and quality control (trim-galore v0.4.1 using `--paired` and `--illumina` options)
- Deduplication (fastuniq v1.1)
- Additional steps for counting reads at each stage.
Run with `snakemake --snakefile Snakereads -j {threads}`

#### Snakedeplete

This Snakefile maps the cleaned reads against the host genome and depletes all the mapping reads
Run with `snakemake --snakefile Snakedeplete -j {threads}`


#### Snakeassemble

This Snakefile runs processed reads through de novo assembly using SPAdes v.3.13
Run with `snakemake --snakefile Snakeassemble -j {threads}`

#### Snakediamond

This Snakefile classifies contigs from the de novo assembly through several steps:
   - BLASTX using Diamond v.2.0.9 against a Diamond2 database including all RefSeq Proteins and their taxonomy information
   - Custom Python scripts to match each contig to the Lowest Common Ancestor of all its BLASTx hits.
   - Custom Python scripts to identify and extract contigs matching a target taxon (e.g. all contigs with arenavirus hits)

Run with `snakemake --snakefile Snakediamond -j {threads}`

Diamond hit results are stored in the Diamond directory
Classified contigs are stored in the Contigs directory
    - Each target taxon has its own subdirectory.


