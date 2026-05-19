test_that("ageTMP_status describes framework scope", {
  expect_match(ageTMP_status(), "age-dependent tumor molecular trajectory")
  expect_match(ageTMP_status(), "CPSA")
})
