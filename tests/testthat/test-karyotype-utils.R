test_that("karyotype parsing round-trips", {
  labels <- c("2.2.1", "3.2.1", "2.3.1")
  mat <- karst:::parse_karyotypes(labels)
  expect_equal(as.character(karst:::format_karyotypes(mat)), labels)
})
