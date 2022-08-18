test_that("is.exstra_db() works", {
  expect_true(is.exstra_db(exstra_known))
  expect_true(is.exstra_db(exstra_wgs_pcr_2))
  expect_false(is.exstra_db(pi))
})

test_that("exstra_db_new_() works", {
  expect_error(exstra_db_new_(5))
  expect_true(inherits(exstra_db_new_(exstra_known$db), "exstra_db"))
  expect_true(inherits(exstra_db_new_(setnames(copy(exstra_known$db), "locus", "disease.symbol")), "exstra_db"))
})

test_that("print.exstra_db() works", {
  expect_output(print(exstra_known), "exstra_db object with 21 loci \\(\\$db\\) of type named")
})

test_that("loci on exstra_db works", {
  expect_true(all(is.element(loci(exstra_known), exstra_known$db$locus)))
  expect_true(all(is.element(exstra_known$db$locus, loci(exstra_known))))
})



# loci_text_info.exstra_d

# loci_normal

# loci_min_exp

# `[.exstra_db`

# copy.exstra_db

# verify.exstra_db()
X_verify <- copy(exstra_known)
setkey(X_verify$db, "motif")
test_that("verify.exstra_db() works", {
  local_edition(2)
  expect_true(verify.exstra_db(exstra_known))
  expect_error(verify.exstra_db(X_verify)) 
  expect_error(verify.exstra_db(6))
})

