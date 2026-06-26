test_that("sceua converges on a two-dimensional sphere", {
  set.seed(1969)
  result <- sceua(
    fn = function(x) sum(x^2),
    lower = c(-5, -5),
    upper = c(5, 5),
    max_evaluations = 5000L,
    kstop = 5L,
    pcento = 1e-8,
    complexes = 5L
  )

  expect_s3_class(result, "sceua")
  expect_length(result$par, 2)
  expect_lt(result$value, 1e-6)
  expect_gt(result$counts, 0)
  expect_gt(result$iterations, 0)
  expect_in(
    result$termination,
    c("objective_convergence", "parameter_convergence", "max_evaluations")
  )
  expect_true(is.data.frame(result$history))
})

test_that("sceua passes extra arguments to the objective", {
  set.seed(1969)
  fn <- function(x, target) sum((x - target)^2)

  result <- sceua(
    fn = fn,
    lower = c(-5, -5),
    upper = c(5, 5),
    target = c(1, 2),
    max_evaluations = 5000L
  )

  expect_lt(sum((result$par - c(1, 2))^2), 1e-2)
})

test_that("sceua validates bound lengths", {
  expect_error(
    sceua(fn = function(x) sum(x^2), lower = c(-5), upper = c(5, 5)),
    "same length"
  )
})

test_that("sceua respects initial point", {
  set.seed(1969)
  result <- sceua(
    fn = function(x) sum(x^2),
    lower = c(-5, -5),
    upper = c(5, 5),
    initial = c(1, 1),
    max_evaluations = 100L
  )

  expect_length(result$par, 2)
  expect_true(result$value < Inf)
})

test_that("sceua inherits R's RNG state", {
  run_sceua <- function() {
    sceua(
      fn = function(x) sum(x^2),
      lower = c(-5, -5),
      upper = c(5, 5),
      max_evaluations = 500L
    )
  }

  # Test that the same seed produces the same result
  set.seed(42)
  first <- run_sceua()
  set.seed(42)
  second <- run_sceua()

  expect_equal(first$par, second$par)
  expect_equal(first$value, second$value)
  expect_equal(first$history, second$history)

  # Test that different seeds produce different results
  set.seed(1234)
  third <- run_sceua()

  expect_false(isTRUE(all.equal(first$par, third$par)))
  expect_false(isTRUE(all.equal(first$value, third$value)))
  expect_false(isTRUE(all.equal(first$history, third$history)))
})
