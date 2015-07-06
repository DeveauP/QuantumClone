# QuantumClone and QuantumCat

R package also available on [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html)
Maintainer: Paul Deveau (paul.deveau at curie.fr)

## Clonal Reconstruction from High-Throughput Sequencing data
Beginners instruction are designed so that you don't see a line of code. Use this if you are not familiar with R programming, and thus that you will use the Graphical User Interface (GUI).

Advanced instructions assume you know the basics of programming (downloading and using packages from CRAN).

This Readme is divided in four parts:

1. Installation instructions (Beginner)

2. Installation instructions (Advanced)

3. Usage (Beginner)

4. Usage (Advanced) 



### Installation instructions (Beginner)

QuantumClone is a method written in R (version 3.1). [R can be downloaded here](http://cran.r-project.org/mirrors.html). 
A user-friendly environment is [RStudio](http://www.rstudio.com/products/rstudio/download/), that can be downloaded to see the data and the code.

In R or RStudio, type (or copy/paste) the following instructions:
> install(ggplot2) 

> install(fpc)

> install(parallel)

> install(doSNOW)

> install(QuantumClone)

You have now successfully installed the QuantumClone package!

It is now time to install and launch the graphical user interface. In R/Rstudio type:
* install(RGtk2)
* library(RGtk2)

Now download the file GUI.R from the master. 
* With R, type: source(/PATH/TO/GUI.R)
* With RStudio, go to Code > Source file and choose GUI.R (*to be added*)

If everything went well, you should see:
![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/GUI.png)
 
*Note*
While the install part can be done once and for all, loading the libraries is mandatory for each session. 

### Installation instructions (Advanced)

The full package is available and is maintained on [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html). 
You can chose to source the GUI but it may affect computation performances.

### Usage (Beginner)
QuantumClone is looking for clones in your samples assuming that there is an evolutionary logic between samples, so you should use data from the same patient for one analysis (either different timepoints, or spatially separated samples, or biological replicates).


QuantumClone requires few informations in the input file:
* The columns in the file MUST be separated by tabulations
* Line 1 should be the column titles (Sample | Chr | Start | Alt | Depth ). An additional argument is required if you do not have a FREEC profile associated to your files: the Genotype. 
* The first column needs to be the name of your sample
* The Chr column contains the chromosome of variant (e.g. "chr2")
* Start is the position of the variant
* Alt is the number of reads supporting the variant
* Depth is the depth of coverage at the position of the variant (number of reads mapped at this position)

While the input file can be as large as you want, the computation time will exponentially grow with the number of variants to be studied. In order to keep computation time reasonable (from a minute to an hour), a reasonable set of mutation is between 100 to 1000 variants.



### Usage (Advanced)
The QuantumClone package is divided in two:
* The clonal reconstruction: QuantumClone / One_step_clustering functions
* The clonal simulation: QuantumCat (not included in the GUI)


