test_that("score_overlap_method()", {
  expect_equal(score_overlap_method("TTCAGGCATT", "CAG"), 6)
  expect_equal(score_overlap_method("TTCAGCAGTT", "CAG"), 6)
  expect_equal(score_overlap_method("TTCAAGCATT", "CAG"), 4)
  expect_equal(score_overlap_method("CAGCAGCAGCA", "CAG"), 11)
})

test_that("score_count_method()", {
  expect_equal(score_count_method("TTCAGGCATT", "CAG"), 2)
  expect_equal(score_count_method("TTCAGCAGTT", "CAG"), 4)
  expect_equal(score_count_method("TTCAAGCATT", "CAG"), 2)
  expect_equal(score_count_method("CAGCAGCAGCA", "CAG"), 9)
})
