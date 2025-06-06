test_that("An example for a real data model fitting", {
data(data)
data(adjacency_matrix)
result=SeasEpi_Par_Est(data,adjacency_matrix,2,2,0.5, 0.5, 1, 0.1, 1, 1, 1, 20, 2,0.2,0.2,5)
expect_type(result, "list")
})
