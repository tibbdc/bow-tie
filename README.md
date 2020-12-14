# bow tie
Determine the bow tie structure by calculating pathways between metabolites using a flux balance analysis (FBA) based approach.
## About

The pipeline was written and tested with Python 3.8. The core libraries essential for the pipeline including: Cobra, Pandas, networkx, and related packages. 

## Installation

1. create bow_tie environment using conda:

```shell
$ conda create -n bow_tie python=3.8
```

2. install related packages using pip:

```shell 
$ conda activate bow_tie
$ pip install cobra
$ pip install networkx
$ pip install openpyxl
$ pip install xlrd
$ pip install ipykernel
$ python -m ipykernel install --user --name bow_tie --display-name "bow_tie"
```

## Steps to reproduce the analysis in the publication

Download all data and analysis code from github (directlt download or use git clone). 

 ```shell
$ cd /file path/project save path/
$ git clone https://github.com/tibbdc/bow-tie.git
```

 All results can be reproduced by executing the Jupyter Python notebooks:

+ bowtie_path.ipynb
  + the main script of determining the bow tie structure by FBA-based approach.

+ bowtie_graph.ipynb
  + the script for determining the bow tie structure by graph

