test_that("data source manifest is well formed", {
  sources <- ageTMP_data_sources(data_dir = tempdir())
  expect_true(all(c("source", "file", "path", "exists") %in% names(sources)))
  expect_true("STable4" %in% sources$source)
  expect_false(any(sources$exists))
})

test_that("sample ID normalization matches paper helper behavior", {
  expect_equal(
    ageTMP_normalize_sample_ids(c("X7316.1052", "A.7316.1099", "P.7316.114")),
    c("7316-1052", "7316-1099", "7316-114")
  )
})
