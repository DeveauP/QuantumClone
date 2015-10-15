[![Build Status](https://travis-ci.org/DeveauP/QuantumClone.svg)](https://travis-ci.org/DeveauP/QuantumClone)
[![CRAN version](http://www.r-pkg.org/badges/version/QuantumClone)](http://www.r-pkg.org/badges/version/QuantumClone)
# QuantumClone and QuantumCat

R package also available on [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html)
Maintainer: Paul Deveau (paul.deveau at curie.fr)

## Clonal Reconstruction from High-Throughput Sequencing data
Beginners instruction are designed so that you don't see a line of code. Use this if you are not familiar with R programming, and thus that you will use the Graphical User Interface (GUI).

Advanced instructions assume you know the basics of programming (downloading and using packages from CRAN).

This Readme is divided in four parts:

[1. Installation instructions (Beginner)](#IIB)

[2. Installation instructions (Advanced)](#IIA)

[3. Usage (Beginner)](#UB)

[4. Usage (Advanced)](#UA) 



### <a name="IIB"></a> Installation instructions (Beginner)

**Preliminary:** This R package has been tested on Windows 8 and Linux. Some issues may appear with OSX setups due to the parallel package apparently. Do not hesitate to contact the maintainer for troubleshooting.

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
* With RStudio, go to Code > Source file and choose GUI.R (**Available for Windows. Linux and OSX GUIs to be added**)

If everything went well, you should see:
![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/GUI.png)
 
 Once you click on the "Validate", the analysis is run and you get a popup notification:
 ![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/popup.png)
 
 When the analysis is complete, an "Exit" button appears, telling you that everything is fine.
 ![alt tag](https://github.com/DeveauP/QuantumClone/blob/master/Images/popup_finished.png)
 
**Note:**
All the libraries are called by the GUI.R file, so there is no need to load them prior to the analysis.



### <a name="IIA"></a>  Installation instructions (Advanced)

The full package is available and is maintained on [CRAN](http://cran.r-project.org/web/packages/QuantumClone/index.html). 
You can chose to source the GUI but it may affect computation performances.

**Note to OSX users:** parallel package seems to be unavailable for R 3.1.2
One option is to download the tarball and use the "R CMD INSTALL QuantumClone.tar.gz" commandline.

### <a name="UB"></a>  Usage (Beginner)
QuantumClone is looking for clones in your samples assuming that there is an evolutionary logic between samples, so you should use data from the same patient for one analysis (either different timepoints, or spatially separated samples, or biological replicates).


QuantumClone requires few informations in the input file:
<a name="SNVD"></a>
* The columns in the file MUST be separated by tabulations
* Line 1 should be the column titles (Sample | Chr | Start | Alt | Depth ). An additional argument is required if you do not have a [FREEC](http://bioinfo-out.curie.fr/projects/freec/) profile associated to your files: the Genotype. 
* The first column needs to be the name of your sample
* The Chr column contains the chromosome of variant (e.g. "chr2")
* Start is the position of the variant
* Alt is the number of reads supporting the variant
* Depth is the depth of coverage at the position of the variant (number of reads mapped at this position)

**Any additional column will not be taken into account for the analysis**

You should have something similar to this:

![alt tag](https://github.com/DeveauP/QuantumClone/raw/master/Images/Example_input.png)


While the input file can be as large as you want, the computation time will exponentially grow with the number of variants to be studied. In order to keep computation time reasonable (from a minute to an hour), a reasonable set of mutation is between **100 to 1000 variants**.

* <a name="FREECD"></a> [FREEC](http://bioinfo-out.curie.fr/projects/freec/) files: list of files corresponding to your samples. It is required if you do not have a Genotype column in your analysis. You should use the "Sample_ratio.txt" file, not the "Sample_ratio_normal.txt" 
* Contamination: fraction of normal cells estimated to contaminate your samples. Needs to be separated by commas (example: 0.1, 0.2)
* Clone range: how many clones should be looked for in the samples? "2:5" means 2 to 5, whereas "2,5" means 2 and 5.
* Save plot: Do you want to save 2D plots?
* Save data: Do you want to keep probabilities and estimated copy numbers in a file?

**If Save plot is ticked you should get two outputs: one in black giving the position of ALL the possibilities for all variants (a), and the second being the result of the clustering for the most likely position of all variants (b).**

(a) ![alt tag](https://github.com/DeveauP/QuantumClone/raw/master/Images/Cellularity1_1.png)
(b) ![alt tag](https://github.com/DeveauP/QuantumClone/raw/master/Images/Cellularity_clustered1_1_1_2.png)

Below is an example of 3D output generated by rgl:
![alt tag](https://github.com/DeveauP/QuantumClone/raw/master/Images/Example_3D.png)

### <a name="UA"></a>  Usage (Advanced)
The QuantumClone package is divided in two:
* [The clonal reconstruction](#CR): QuantumClone / One_step_clustering functions
* [The clonal simulation](#CS): QuantumCat (not included in the GUI)

#### <a name="CR"></a> Clonal reconstruction
One_step_clustering() has several parameters required (some have default configuration):
> One_step_clustering(SNV_list, FREEC_list = NULL, contamination,
  nclone_range = 2:5, clone_priors = NULL, prior_weight = NULL,
  maxit = 1 , preclustering = T, simulated = F, epsilon = 5 * (10^(-3)),
  save_plot = T, ncores = 1, plot_3D = F, plot_3D_before_clustering = F,
  restrict.to.AB = F, output_directory = NULL)

* SNV_list: list of dataframes. See [previous section](#SNVD) for description.
* FREEC_list: list of outputs from FREEC (in the same order as the SNV list). See [here](#FREECD) for added information.
* contamination: Numeric vector giving the fraction of normal cells in each sample. Is linked to the cellularity by contamination = 1 - Cellularity
* nclone_range: number of clones to look for in the samples
* clone_priors: list of vectors giving the position of the clones in each samples (if know from previous analysis)
* prior_weight : fraction of variants belonging to a clone (if known from previous analysis)
* maxit : number of iterations to run per condition. The output will take the maximal maximum likelihood on all iterations.
* preclustering : should kmeans (fpc package) be used to give priors to the alorithm?
* simulated : is the data generated by QuantumCat? It does not change the parameters, but will attribute shapes to different chromosomes in the plots. (see [QuantumCat](#CS) for more information)
* epsilon : stop condition for the EM.
* save_plot : save the 2D plots in a folder with the patient name/output_directory.
* ncores: number of CPUs on which to distribute calculations (used if high number of variants)
* restrict.to.AB : should the clustering be done only on AB regions?
* output_directory : directory in which the plots will be saved (if NULL, will create a directory with the patient name)


#### <a name="CS"></a> Clonal simulation
This part is about generating data to test clonal reconstruction algorithms. Its core is the QuantumCat function. It will generate data for a single cancer that can be sequenced multiple times (either spatially separated or different timepoints). It thus assumes that there is an evolutionary history between samples. The "Chr" columns stores the information of the clonal attribution.
> QuantumCat(number_of_clones, number_of_mutations, ploidy = 2, depth = 100,
  number_of_samples = 2, Random_clones = F, contamination = NULL)

* number_of_clones : How many clones should exist in total. For example, 5 clones in 2 samples can be distributed in the following way: 1 specific of sample 1, 1 specific of sample 2 and 3 shared between sample 1 and 2.
* number_of_mutations : How many variants should be used for the clustering. Some algorithms reported that an increase in the number of variants decreased the clustering quality, which does not seem to be the case here. It affects the computing time however.
* ploidy : if numeric, it will generate a poisson distribution with mean the ploidy. Accepted inputs can be "disomic", "AB", "AAB", "A", etc.
* depth : what is the sequencing depth? Depth of a variant will be generated according to a negative binomial distribution, which characteristics have been generated by fitting to our data from whole genome sequencing.
* number_of_samples = How many samples should be generated?
* Random_clones: if the number of clones should be generated randomly (sampled from 2:5)
* contamination: estimation of the contamination by normal cells

For multiple testings, and calculation of the Normalized Mutual Information (NMI), see Multitest() and statistics_on_Multitest()

###Perspectives
The doSnow library for parallel computing has been replaced by doParallel in the latest R version. This is known and a new version will soon be available taking into account this change.

### Acknowledgments
Many thanks to the contributors of this work: my supervisors, Elodie for the features improvement and Linux debugging, Matahi for the OSX feedback, and more generally to the U830 & U900 people. This work had been funded by the Ministere de l'Enseignement Supérieur de la Recherche (AMX grant).
