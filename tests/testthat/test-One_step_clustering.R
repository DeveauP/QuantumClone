context("Clustering quality and options")

test_that("Check optim",{
  expect_identical(Compute_NMI(
    QuantumClone::One_step_clustering(Input_Example,contamination = c(0,0),nclone_range = 3:5,optim = "default",save_plot= FALSE,preclustering = "FLASH")
  ), 1)
  
})
test_that("Check optimx",{
  expect_identical(Compute_NMI(
    QuantumClone::One_step_clustering(Input_Example,contamination = c(0,0),
                                      nclone_range = 3:5,
                                      optim = "optimx",
                                      save_plot = FALSE,
                                      preclustering = "FLASH")
  ), 1)
  
})

test_that("Check compound, diploid",{
  QC<-QuantumClone(Input_Example,contamination = c(0,0), nclone_range = 3:5,
                   save_plot = FALSE, optim = "compound",
                   Init =1 ,epsilon = 0.01,
                   preclustering = "FLASH",
                   keep.all.model = FALSE,
                   model.selection = 2)
   expect_identical(Compute_NMI(QC),1)
})