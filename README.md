# YaMaS (YOLO lab Microbiome System)

YaMaS is a package designed to easily download DNA datasets from the NCBI SRA website. It is developed by the YOLO lab team, and is designed to be simple, efficient, and easy to use for non-programmer users.

## Installation

To install YaMaS, you can use pip:

```
pip install yamas
```

## Dependencies
Before proceeding with the installation of YaMaS, please make sure that all dependencies are fulfilled. In case any of the dependencies are missing, the program will not run as expected. Please refer to the installation instructions and ensure that all requirements are met before proceeding.
- YaMaS should be downloaded in a [qiime2](https://docs.qiime2.org/2023.2/) enviorment.
- [SRA-toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) package should be downloaded in the enviorment.
- [Entrez](https://anaconda.org/bioconda/entrez-direct) package should be downloaded in the enviorment.
- Exporting a project requires a downloaded [classifier file](https://data.qiime2.org/<qiime2-version>/common/gg-13-8-99-nb-classifier.qza).

## Getting Started

YaMaS provides an easy-to-use interface in the terminal.

To download the visualization file for multiple projects and save each file to a separate folder, use the following command:
```
yamas --download PRJEB12345
```
To export an OTU (Operational Taxonomic Unit), taxonomy, and phylogeny tree for a single project, use the following command:
```
yamas --export project_path 0 40 classifier_file otu_output_file taxonomy_output_file 16
```
In the --export command, the project_path argument should be the path to the project directory, classifier_file should be the path to the trained classifier file, and otu_output_file and taxonomy_output_file should be the paths to the output files. The 16 argument specifies the number of threads to use for parallel processing, which can speed up the export process (default is 12).
