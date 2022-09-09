test_that("add_alpha() works", {
  expect_equal(add_alpha_("blue", 0.5), c(blue = "#0000FF80"))
  expect_equal(add_alpha_("#0000FF", 0.5), c("#0000FF" = "#0000FF80"))
  expect_equal(add_alpha_("hotpink", 0.2), c(hotpink ="#FF69B433"))
  expect_equal(add_alpha_("#FF69B4", 0.2), c("#FF69B4" = "#FF69B433"))
  expect_equal(add_alpha_(c(), 0.2), c())
  expect_error(add_alpha_("not a colour"))
  expect_error(add_alpha_(), "Please provide a vector of colours")
})
