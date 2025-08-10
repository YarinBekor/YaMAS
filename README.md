# YaMAS (YOLO Microbiome Analysis System)

### Table of Contents
- [Pipeline Overview](#-pipeline-overview)
- [Installation & Dependencies](#-installation--dependencies)
- [Flags & Options](#-flags--options)
- [Output Structure](#-output-structure)
- [Downloading a project](#downloading-a-project)
    - [Download from NCBI SRA](#download-from-ncbi-sra)
    - [Continue data downloading](#continue-data-downloading)
    - [Download from ENA](#download-from-ena)
    - [Download using fastq files](#download-using-fastq-files)
- [Exporting a project (only for 16S/18S)](#exporting-a-project-only-for-16s18s)
- [Arguments and configurations](#arguments-and-configurations)


YaMAS is a package designed to easily download DNA datasets from the NCBI SRA,ENA and qiita websites. It is developed by the [YOLO lab team](https://yolo.math.biu.ac.il), and is designed to be simple, efficient, and easy to use for non-programmers.

## Pipeline Overview

The YaMAS pipeline consists of several stages, depending on the sequencing type:

**Supported input sources:**  
- **Project ID** from SRA/ENA/Qiita (automatic download)  
- **Local FASTQ file** (process without download)  
- **Existing FASTQ folder** (`--continue_from_fastq`)  
- **Existing project folder** (`--continue_from`)  

**For Shotgun datasets:**  
1. **Download** dataset from SRA/ENA/Qiita  
2. **Preprocessing** (soon - optional host removal, quality control)  
3. **MetaPhlAn** – Taxonomic profiling  
4. **HUMAnN** *(if `--pathways` flag is set)* – Functional profiling and pathway analysis  
5. **Export & Visualization** – Generation of merged abundance tables and plots  

**For 16S/18S datasets:**  
1. **Download** dataset  
2. **QIIME2 processing** – Denoising, taxonomic classification  
3. **Export & Visualization**  

> **Note:** HUMAnN integration is available **only** for Shotgun datasets and runs immediately after MetaPhlAn.

##  Installation & Dependencies
Before proceeding with the installation and execution of YaMAS, please ensure that you have a clean environment set up on your system, with all dependencies installed. To create one, follow the steps below:

1. Create a new [qiime2](https://docs.qiime2.org/2023.2/install/native/) environment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html). Make sure you name it 'qiime2'.
2. Download the [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) and [Entrez](http://bioconda.github.io/recipes/entrez-direct/README.html) packages to the environment.
3. Download the [metaphlan](https://github.com/biobakery/biobakery/wiki/metaphlan4) package. Make sure the database works properly before proceeding.
4. Exporting a 16S project requires a downloaded [classifier file](https://data.qiime2.org/2022.8/common/gg-13-8-99-nb-classifier.qza).
5. Get YaMAS [ready](https://github.com/YarinBekor/YaMAS#get-yamas-ready). 

### Step-by-Step: Setting Up the YaMAS Environment

1. **Install QIIME2 environment**  

```
wget https://data.qiime2.org/distro/core/qiime2-2023.2-py38-linux-conda.yml
conda env create -n qiime2-2023.2 --file qiime2-2023.2-py38-linux-conda.yml
conda activate qiime2-2023.2
```

2. **Install required tools**  
```
conda install -c conda-forge mamba
```
Install SRA Tools via Mamba
```
mamba install -c bioconda sra-tools
#Optional check:
which prefetch
which fasterq-dump
prefetch --version
fasterq-dump --version
```
Install Entrez Direct
```
mamba install -c bioconda entrez-direct
#Optional check:
which esearch
```
Install MetaPhlAn (v3.0.7)
```
mamba install -c bioconda metaphlan=3.0.7
#Optional check:
metaphlan --version
```
Install HUMAnN 
```
mamba install -c biobakery humann
```

3. **Download & configure the MetaPhlAn database**
```
metaphlan --install --bowtie2db /path/to/db --index mpa_v30_CHOCOPhlAn_201901
export METAPHLAN_BOWTIE2_DB=/path/to/db
echo 'export METAPHLAN_BOWTIE2_DB=/path/to/db' >> ~/.bashrc
```

4. **Set HUMAnN database**
```
humann_config --update database_folders nucleotide /path/to/db/chocophlan

humann_config --update database_folders protein /path/to/db/uniref

humann_config --update database_folders utility_mapping /path/to/db/utility_mapping
```

5. **Install YaMAS**
```
pip install YMS
```

Get YaMAS ready:
```
yamas --ready <operating_system_type> 
```
Arguments:
- operating_system_type: Ubuntu/CentOS

> Pay attention to the output of the command.    
If the environment is ready, you will need to run one more command.    
If not, follow the output guidelines.  

**You’re all set!**  

## Flags & Options
| Flag | Description |
|------|-------------|
| `--download <PROJECT_ID>` | Download a dataset from SRA/ENA/Qiita |
| `--type <16S/18S/Shotgun>` | Type of sequencing data |
| `--as_single` | Treat paired-end reads as single-end |
| `--pathways yes/no` | Enable HUMAnN for pathway profiling (Shotgun only) |
| `--threads <N>` | Number of threads to use |
| `--continue_from_fastq <ID> <PATH> <TYPE>` | Continue processing from an existing FASTQ folder |
| `--continue_from <ID> <PATH> <TYPE>` | Continue processing from an existing dataset folder |

---

## Output Structure
After running YaMAS, the project folder will contain:

```
<PROJECT_ID>-<DATE>_<TIME>/
│
├── sra/                # Raw SRA files
├── fastq/              # FASTQ files 
├── qza/                # QIIME2 artifacts (16S/18S)
├── vis/                # Visualization files
├── export/             # Exported tables and merged results
└── humann_results/     # (If --pathways is set) HUMAnN output files
```
**HUMAnN outputs include:**  
- `*_pathabundance.tsv` – Normalized pathway abundance per sample  
- `*_pathcoverage.tsv` – Pathway coverage per sample  
- `*_pathabundance_stratified.tsv` – Stratified pathway abundance by species  

---

## Downloading a project
To download a project from **NCBI SRA** or from **ENA, qiita**, use the one of the following templates:    

### Download from NCBI SRA
```
yamas --download <dataset_id> --type <data_type> --pathways <pathways>
```
Arguments:
- dataset_id: the dataset id from the NCBI SRA website. For example: PRJEB01234
- data_type: choose one of the following types: 16S / 18S / Shotgun
- pathways: Generate HUMAnN pathways tables. choose: yes / no 

### Continue data downloading  
1. Continue downloading project **after** downloading SRA **before** converting to .fastq.    
Use the following command:
```
yamas --continue_from_fastq <dataset_id> <project_path> <data_type> --pathways <pathways>
```
Arguments:
- dataset_id: the dataset id from the NCBI SRA website. For example: PRJEB01234
- project_path: path to the project directory (created by YaMAS, if you started downloading data in the past).
- data_type: choose one of the following types: 16S / 18S / Shotgun
- pathways: Generate HUMAnN pathways tables. choose: yes / no
    

2. Continue downloading project **after** downloading SRA and **after** converting them to .fastq.  
Use the following command:
```
yamas --continue_from <dataset_id> <project_path> <data_type> --pathways <pathways>
```
Arguments:
- dataset_id: the dataset id from the NCBI SRA website. For example: PRJEB01234
- project_path: path to the project directory (created by YaMAS, if you started downloading data in the past).
- data_type: choose one of the following types: 16S / 18S / Shotgun
- pathways: Generate HUMAnN pathways tables. choose: yes / no


## Download from ENA
```
yamas --qiita <preprocessed_fastq_path> <metadata_path> <data_type>
```
Arguments:
All can be found in https://qiita.ucsd.edu/   
- data_type : choose one of the following types: 16S / 18S 
- Where preprocessed fastq can be found?    
    Click the study description --> in the graph click on 'demultiplexed' --> scroll down and download 'preprocessed fastq' --> rename the file to be: "forward.fastq.gz"
- Where metadata can be found?
    Click the study description --> download 'Prep info' --> rename the file to be: "metadata.tsv"
- The new data will be created in the folder of the fastq and metadata, so it is recommended to be organized.

### Download using fastq files
```
yamas --fastq <preprocessed_fastq_path> <barcode_path> <metadata_path> <data_type>
```
Arguments:
- preprocessed_fastq_path: path to the preprocessed fastq file. rename the file to be: "preprocessed_fastq_path"
- barcode_path: path to the barcode file. rename the file to be: "barcodes.fastq.gz"
- metadata_path: path to the metadata file. rename the file to be: "metadata.tsv". The metadata should contains column names: "barcode".
- data_type: choose one of the following types: 16S / 18S / Shotgun

---

## Exporting a project (only for 16S/18S)
To export an OTU (Operational Taxonomic Unit), taxonomy, phylogeny tree and a tree.nwk for a single project, use the following command:
```
yamas --export <project_path> <data_type> <start> <end> <classifier_file> <threads>
```
Arguments:
- project_path: path to the project directory (created by YaMAS in the previous step).
- data_type: choose one of the following types: 16S / 18S / Shotgun
- classifier_file: path to the trained classifier file. 
- start & end: choose graph edges. 
- threads: specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).


## Arguments and configurations
1. config: You can add a configuration file in order to save the data in a different folder, and change other configurations. 
2. verbose: To get more information about a downloading process, use the verbose option (this is highly recommended).
3. Listing more than one project will download them one by one into different folders.


## Cite us
If you are using our package, YaMAS for **any** purpose, please cite us; Shtossel Oshrit, Sondra Turjeman, Alona Riumin, Michael R. Goldberg, Arnon Elizur, Yarin Bekor, Hadar Mor, Omry Koren, and Yoram Louzoun. "Recipient-independent, high-accuracy FMT-response prediction and optimization in mice and humans." Microbiome 11, no. 1 (2023): 181. https://link.springer.com/article/10.1186/s40168-023-01623-w