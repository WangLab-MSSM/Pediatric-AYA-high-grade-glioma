test_that("normal reference loader fails clearly when package data are absent", {
  expect_error(
    ageTMP_load_normal_reference(path = file.path(tempdir(), "missing.rds")),
    "Normal-reference data were not found"
  )
})
