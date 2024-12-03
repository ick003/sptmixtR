test_that("generate data works", {
  x <- generate_sptmixtdata()
  expect_type(x, "list")
})
