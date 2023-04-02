# YaMAS (YOLO Microbiome Analysis System)

YaMAS is a package designed to easily download DNA datasets from the NCBI SRA website. It is developed by the [YOLO lab team](https://yolo.math.biu.ac.il), and is designed to be simple, efficient, and easy to use for non-programmers.

## Installation

To install YaMAS, you can use pip:

```
pip install yamas
```

## Dependencies
Before proceeding with the installation and execution of YaMAS, please ensure that you have a clean environment set up on your system, with all dependencies installed. To create one, follow the steps below:
1. create a new [qiime2](https://docs.qiime2.org/2023.2/install/native/) environment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html). Make sure you name it 'qiime2'.
2. download the [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) and [Entrez](http://bioconda.github.io/recipes/entrez-direct/README.html) packages to the environment.
3. Exporting a project requires a downloaded [classifier file](https://data.qiime2.org/2022.8/common/gg-13-8-99-nb-classifier.qza).

You are now ready to run and install YaMAS in the newly created and activated qiime2 environment.
## Getting Started

YaMAS provides an easy-to-use interface in the terminal.

To download a project, use the following command:
```
yamas --download PRJEB01234
```
Listing more than one project will download them one by one into different folders.

To export an OTU (Operational Taxonomic Unit), taxonomy, and phylogeny tree for a single project, use the following command:
```
yamas --export <project_path> <trim> <trunc> <classifier_file> <threads>
```
Arguments:
- project_path: path to the project directory (created by YaMAS in the previous step).
- classifier_file: path to the trained classifier file. 
- trim & trunc: choose graph edges. 
- threads: specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).

