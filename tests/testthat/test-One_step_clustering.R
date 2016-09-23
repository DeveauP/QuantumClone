context("Clustering quality and options")

test_that("Check optim",{
  expect_identical(Compute_NMI(
    One_step_clustering(Input_Example,contamination = c(0,0),nclone_range = 4,optim = "default",save_plot= FALSE)
  ), 1)
  
})
test_that("Check optimx",{
  expect_identical(Compute_NMI(
    One_step_clustering(Input_Example,contamination = c(0,0),nclone_range = 4,optim = "optimx",save_plot = FALSE)
  ), 1)
  
})

