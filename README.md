# bow tie
Determine the bow tie structure by calculating pathways between metabolites using a flux balance analysis (FBA) based approach.
## About

The pipeline was written and tested with Python 3.8. The core libraries essential for the pipeline including: Cobra, Pandas, and related packages. 

## Software

The packages used to run the code in the pipeline was listed in requirements.txt. To install the requirements using pip, run the following at command-line:

```shell
$ pip install -r requirements.txt
```

To create a stand-alone environment named bow_tie with Python 3.8 and all the reqiured package versions, run the following:

```shell
$ conda create -n bow_tie python=3.8 
$ source activate bow_tie
$ pip install -r requirements.txt
```

You can read more about using conda environments in the [Managing Environments](http://conda.pydata.org/docs/using/envs.html) section of the conda documentation. 

## Steps to reproduce the analysis in the publication

 All results can be reproduced by executing the Jupyter Python notebooks:

+ bowtie_path.ipynb
  + the main script of determining the bow tie structure by FBA-based approach.

+ bowtie_graph.ipynb
  + the script for determining the bow tie structure by graph

