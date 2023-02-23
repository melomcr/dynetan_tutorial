# Dynamic Network Analysis - Tutorial
A series of tutorials for the application of Dynamical Network Analysis using Generalized Correlations.

The python module repository can be found at https://github.com/melomcr/dynetan

A full module documentation can be found at https://dynamical-network-analysis.readthedocs.io/en/latest/index.html

This tutorial was prepared and released along with the publication:

* [Generalized correlation-based dynamical network analysis: a new high-performance approach for identifying allosteric communications in molecular dynamics trajectories. J. Chem. Phys. (2020). DOI: 10.1063/5.0018980](https://doi.org/10.1063/5.0018980)

## Tutorial index

**Long-form tutorial introducing the method and showcasing extensive analysis**

This tutorial will introduce the researcher to the fundamentals of Dynamical Network
Analysis and to the entire range of functionalities provided by this python module.
The system used as an example is a single protein enzyme bound to its substrate, 
and its simulations are used as a starting point for correlation calculation, 
community analysis, optimal path determination, comparisons between cartesian and 
network distances, and finally interactive visualization.
The first file is a python notebook that explains the method as it processes the 
trajectory. The second file is a jupyter notebook that exemplifies a pipeline for 
interactive data analysis.
The tutorial concludes with the creation of files for VMD visualization and 
rendering of publication-quality images.

* Tutorial-Step_1-ProcessTrajectory.ipynb
* Tutorial-Step_2-AnalysisAndPlots.ipynb

**Command-Line-Interface version of long-form tutorial**

These files are adaptations of the notebooks in the tutorial above which were 
prepared to facilitate the remote execution of trajectory processing and analysis.

* Tutorial-Command-Line-Interface-Step_1.py
* Tutorial-Command-Line-Interface-Step_2.py

**Single-protein analysis**

This tutorial will cover the common case of a single protein being simulated for
community analysis and network visualization.
The first file is a python script that can be adapted to run in a remote cluster
through a command line interface (CLI). The second file is a jupyter notebook that
serves as a starting point for interactive data analysis. The tutorial concludes
with the creation of files for VMD visualization and rendering of publication-quality 
images.

* Tutorial-Single-Protein-CLI-Step_1.py
* Tutorial-Single-Protein-Step_2.ipynb

**Non-canonical residue and ligands**

This tutorial shows how a researcher can apply the Dynamical Network Analysis 
technique to a system containing non-canonical protein residues, lipids,
carbohydrates, and ligands such as drugs. It shows how interactive visualizations 
and specialized module functions can be used to prepare a network representation 
of complex ligands. 

* Tutorial-Non-Canonical-and-Non-Proteic-Residues.ipynb

## Installing Requirements

There are two main ways to install all dependencies to execute this tutorial: The first is by creating a Conda environment, the second is by manually installing system packages, and Python and R packages. We provide examples for both method below.

Building a Conda Environment
----------------

The easiest way to install all requirements for the execution of this tutorial is to create a Conda environment. Inside the folder *CondaEnv* you will find a recipe for an environment with all required Python and R packages. 
To create the environment using a terminal (a command line interface), simply go to the *CondaEnv* folder and run the following command:

    conda env create -f environment.yml

Custom Installation
----------------

**System Packages**

For this tutorial, you will need an R installation (usually provided by system packages `R-core` and `R-core-devel`) and system packages for the "gdata" package (usually provided by system packages `perl-Compress-Raw-Zlib` and `perl-Spreadsheet-XLSX`). In Fedora 32, all can be installed with the following command:

    dnf install perl-Compress-Raw-Zlib perl-Spreadsheet-XLSX R-core-devel

**Python Packages**

To run the tutorial Jupyter Notebooks, you will need the `dynetan` package and `rpy2`. The following commands will install all required python packages, assuming you have met all system requirements.

    pip install dynetan
    pip install rpy2

To activate notebook widgets (from the ipywidgets package) and visualization (from the nglview package), you may need to execute the following commands:

    jupyter nbextension enable --py widgetsnbextension
    jupyter-nbextension enable nglview --py --user

**R Packages**

To install R packages, execute the following commands in your R terminal:

    list.of.packages <- c("data.table", "ggplot2", "ggrepel", "gdata", "RColorBrewer", "colorRamps", "rPref")

    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if(length(new.packages)) install.packages(new.packages, repos='http://cran.wustl.edu/', dependencies = TRUE)

