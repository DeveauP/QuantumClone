We are committed to reproducible research. To achieve it, we decided to publish all tests scripts used to evaluate QuantumClone and compare it to other algorithms.
QuantumCat was written with this in mind, and we would like to make it easier for everyone to reproduce our results.

# Comparing different conditions
The first folder "Rscripts" contains scripts that can be used to evaluate the clustering quality and time in different settings.

# Comparing QuantumClone to other algorithms
We decided to compare QuantumClone to sciClone [(version 1.0.7)](https://github.com/genome/sciclone), pyClone [(version 0.12.9)](https://bitbucket.org/aroth85/pyclone/wiki/Home), and to a k-medoid clustering of variant allele frequencies [(pamk)](https://cran.r-project.org/web/packages/fpc/fpc.pdf).
All tests show the comparison in clustering quality defined by the Normalized Mutual Information [(NMI)](http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html), and also give an estimation for the computing time.
All settings are tested by 50 randomly generated tumors (which can contain multiple sample each), and for this reason the computational time can become high (a few days worth of testing).

All suggestions are welcomed to improve testing and algorithms.