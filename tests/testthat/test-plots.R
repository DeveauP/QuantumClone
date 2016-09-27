context("Plots")

test_that("plot_QC_out returns a ggplot object", {
  g<-QuantumClone::plot_QC_out(QC_output,Sample_names = c("Diag","Relapse"))
  expect_identical(class(g),c("gg","ggplot"))
}
)

test_that("plot_QC_out returns a ggplot object", {
  g<-QuantumClone::evolution_plot(QC_output,Sample_names = c("Diag","Relapse"))
  expect_identical(class(g),c("gg","ggplot"))
}
)