test_that("An example for a simulation study", {
  result=SeasEpi_Sim_Par_Est(5,5,10,30,0.7, 0.7, -1, 0.1, 0,40, 50,0.6, 5, 5, 10, 3,0.2,0.2,5)
  expect_type(result, "list")
})
