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

**Preliminary:** This R package has been tested on Windows 8 and Linux. Some issues may appear with Mac setups due to the rgl package apparently. Do not hesitate to contact the maintainer for troubleshooting.

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

Now download the file GUI.R from the master. 
* With R, type: source(/PATH/TO/GUI.R)
* With RStudio, go to Code > Source file and choose GUI.R (**to be added**)

If everything went well, you should see:
![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/GUI.png)
 
**Note:**
All the libraries are called by the GUI.R file, so there is no need to load them prior to the analysis.

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

**Any additional column will not be taken into account for the analysis**

You should have something similar to this:
![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/Example_input.png)


While the input file can be as large as you want, the computation time will exponentially grow with the number of variants to be studied. In order to keep computation time reasonable (from a minute to an hour), a reasonable set of mutation is between **100 to 1000 variants**.

* FREEC files: list of files corresponding to your samples. It is required if you do not have a Genotype column in your analysis.
* Contamination: fraction of normal cells estimated to contaminate you samples. Needs to be separated by commas (example: 0.1, 0.2)
* Clone range: how many clones should be looked for in the samples? "2:5" means 2 to 5, whereas "2,5" means 2 and 5.
* Save plot: Do you want to save 2D plots?
* Save data: Do you want to keep probabilities and estimated copy numbers in a file?

**If Save plot is ticked you should get two outputs: one in black giving the position of ALL the possibilities for all variants (a), and the second being the result of the clustering for the most likely position of all variants (b).**

(a) ![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/Cellularity1_1.png)
(b) ![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/Cellularity_clustered1_1_1_2.png)



### Usage (Advanced)
The QuantumClone package is divided in two:
* The clonal reconstruction: QuantumClone / One_step_clustering functions
* The clonal simulation: QuantumCat (not included in the GUI)


