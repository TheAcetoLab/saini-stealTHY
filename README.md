# StealTHY : CRISPR-Cas9 screen

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr



## Data pre-processing

Raw data is available at Gene Expression Omnibus (GEO, NCBI; accession numbers [GSE268044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268044),  [GSE268045](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268045) and [GSE268588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268588) ). Data processing is computationally expensive and is not covered in this repository. We provide description of the data pre-processing workflow together with software version in the original publication. Processed data and large result files are  archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11374563.svg)](https://doi.org/10.5281/zenodo.11374563).



##  Data and code availability

To reproduce our analysis, first clone source code from the [GitHub repository](https://github.com/TheAcetoLab/saini-stealTHY). This repository is also archived at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11288553.svg)](https://doi.org/10.5281/11288553)

    git clone https://github.com/TheAcetoLab/saini-stealTHY

Next, download processed data and output data deposited in [Zenodo](https://doi.org/10.5281/zenodo.11374563) into the cloned project folder and untar the files.

    for file in *.tar.gz; do tar xzvf "${file}" && rm "${file}"; done
