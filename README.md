# QuantumClone and QuantumCat

R package also available on [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html)
Maintainer: Paul Deveau (paul.deveau at curie.fr)

## Clonal Reconstruction from High-Throughput Sequencing data
Beginners instruction are designed so that you don't see a line of code. Use this if you are not familiar with R programming.
Advanced instructions assume you know the basics of programming (downloading and using packages from CRAN)

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
* With RStudio, go to Code > Source file and choose GUI.R

If everything went well, you should see:
![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/GUI.png)
 
*Note*
While the install part can be done once and for all, loading the libraries is mandatory for each session. 

### Installation instructions (Advanced)

The full package is available and is maintained on: [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html). 

### Usage
