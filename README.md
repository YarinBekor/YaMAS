# YaMAS (YOLO Microbiome Analysis System)

YaMAS is a package designed to easily download DNA datasets from the NCBI SRA and ENA website. It is developed by the [YOLO lab team](https://yolo.math.biu.ac.il), and is designed to be simple, efficient, and easy to use for non-programmers.

## Dependencies
Before proceeding with the installation and execution of YaMAS, please ensure that you have a clean environment set up on your system, with all dependencies installed. To create one, follow the steps below:
1. Create a new [qiime2](https://docs.qiime2.org/2023.2/install/native/) environment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html). Make sure you name it 'qiime2'.
2. Download the [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) and [Entrez](http://bioconda.github.io/recipes/entrez-direct/README.html) packages to the environment.
3. Download the [metaphlan](https://github.com/biobakery/biobakery/wiki/metaphlan4) package. Make sure the database works properly before proceeding.
4. Exporting a 16S project requires a downloaded [classifier file](https://data.qiime2.org/2022.8/common/gg-13-8-99-nb-classifier.qza).
5. Get YaMAS [ready](https://github.com/YarinBekor/YaMAS#get-yamas-ready). 

You are now ready to run and install YaMAS in the newly created and activated qiime2 environment.
## Installation

To install YaMAS, you can use pip:

```
pip install YMS
```

## Getting Started- NCBI SRA

YaMAS provides an easy-to-use interface in the terminal.    
To download a project from **NCBI SRA**, use the one of the following templates:    

### Get YaMAS ready
```
yamas --ready <operating_system_type> 
```
Arguments:
- operating_system_type: Ubuntu/CentOS

Pay attention to the output of the command.    
If the environment is ready, you will need to run one more command.    
If not, follow the output guidelines.    

### 16S/18S dataset
```
yamas --download PRJEB01234 --type 16S/18S 
```
To export an OTU (Operational Taxonomic Unit), taxonomy, and phylogeny tree for a single project, use the following command:
```
yamas --export <project_path> <data_type> <start> <end> <classifier_file> <threads>
```
Arguments:
- project_path: path to the project directory (created by YaMAS in the previous step).
- data_type: choose one of the following types: 16S / 18S 
- classifier_file: path to the trained classifier file. 
- start & end: choose graph edges. 
- threads: specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).

### Shotgun dataset
```
yamas --download PRJEB01234 --type Shotgun 
```

### Continue data downloading  
Continue downloading project **after** downloading SRA **before** converting to .fastq.    
Use the following command:
```
yamas --continue_from_fastq <dataset_id> <project_path> <data_type>
```
Arguments:
- project_path: path to the project directory (created by YaMAS, if you started downloading data in the past).
- data_type: choose one of the following types: 16S / 18S / Shotgun    
    

Continue downloading project **after** downloading SRA **and** converting them to .fastq.  
Use the following command:
```
yamas --continue_from <dataset_id> <project_path> <data_type>
```
Arguments:
- project_path: path to the project directory (created by YaMAS, if you started downloading data in the past).
- data_type: choose one of the following types: 16S / 18S / Shotgun

## Getting Started- ENA

YaMAS provides an easy-to-use interface in the terminal.    
To download a project from **ENA,qiita**, use the one of the following templates:    

### Get YaMAS ready
```
yamas --ready <operating_system_type> 
```
Arguments:
- operating_system_type: Ubuntu/CentOS

Pay attention to the output of the command.    
If the environment is ready, you will need to run one more command.    
If not, follow the output guidelines.    

### 16S/18S dataset
```
yamas --qiita <preprocessed_fastq_path> <metadata_path> <data_type>
```
Arguments:
All can be found in https://qiita.ucsd.edu/    
- Where preprocessed fastq can be found?    
    Click the study description --> in the graph click on 'demultiplexed' --> scroll down and download 'preprocessed fastq' --> rename the file to be: "forward.fastq.gz"
- Where metadata can be found?
    Click the study description --> download 'Prep info' --> rename the file to be: "metadata.tsv"
- The new data will be created in the folder of the fastq and metadata, so it is recommended to be organized.

To export an OTU (Operational Taxonomic Unit), taxonomy, and phylogeny tree for a single project, use the following command:
```
yamas --export <project_path> <data_type> <start> <end> <classifier_file> <threads>
```
Arguments:
- project_path: path to the project directory (created by YaMAS in the previous step).
- data_type: choose one of the following types: 16S / 18S 
- classifier_file: path to the trained classifier file. 
- start & end: choose graph edges. 
- threads: specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).


## Arguments and configurations
1. config: You can add a configuration file in order to save the data in a different folder, and change other configurations. 
2. verbose: To get more information about a downloading process, use the verbose option (this is highly recommended).
3. Listing more than one project will download them one by one into different folders.


## Cite us
If you are using our package, YaMAS for **any** purpose, please cite us; Shtossel, Oshrit, Sondra Turjeman, Alona Riumin, Michael R. Goldberg, Arnon Elizur, Yarin Bekor, Hadar Mor, Omry Koren, and Yoram Louzoun. "Recipient-independent, high-accuracy FMT-response prediction and optimization in mice and humans." Microbiome 11, no. 1 (2023): 181. https://link.springer.com/article/10.1186/s40168-023-01623-w