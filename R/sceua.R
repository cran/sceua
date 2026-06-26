#' Minimize a function with SCE-UA
#'
#' @description
#' Find the parameter set that minimizes an objective function using the
#' Shuffled Complex Evolution - University of Arizona (SCE-UA) algorithm
#' (Duan et al., 1992).
#'
#' @param fn Function to minimize. Must accept a single numeric vector of
#'   parameters and return a scalar numeric value.
#' @param lower Numeric vector of lower bounds. Must have the same length as
#'   `upper`.
#' @param upper Numeric vector of upper bounds. Must have the same length as
#'   `lower`.
#' @param initial Optional initial parameter vector. If provided, it is
#'   included in the initial population.
#' @param max_evaluations Maximum number of function evaluations.
#' @param kstop Number of shuffling loops over which the objective value must
#'   change by `pcento` before convergence.
#' @param pcento Objective convergence threshold.
#' @param complexes Number of complexes in the initial population.
#' @param points_per_complex Number of points in each complex. Defaults to
#'   `2 * n + 1` where `n` is the number of parameters.
#' @param simplex_size Number of points in each sub-complex. Defaults to
#'   `n + 1`.
#' @param evolution_steps Number of evolution steps allowed for each complex
#'   before shuffling. Defaults to `points_per_complex`.
#' @param min_complexes Minimum number of complexes required. Defaults to
#'   `complexes`.
#' @param parameter_epsilon Parameter convergence threshold.
#' @param ... Additional arguments passed to `fn`.
#'
#' @details
#' The R wrapper draws the internal SCE-UA seed from R's global random number
#' generator. Call `set.seed()` before `sceua()` for reproducible results.
#'
#' @returns
#' An object of class `sceua`: a list with components:
#' - `par`: best parameter vector.
#' - `value`: objective value at `par`.
#' - `counts`: number of function evaluations.
#' - `iterations`: number of shuffling loops.
#' - `termination`: reason for termination.
#' - `history`: a `data.frame` with one row per shuffling loop.
#'
#' @references
#' Duan, Q., Sorooshian, S., and Gupta, V.K., 1992. Effective and efficient
#' global optimization for conceptual rainfall-runoff models.
#' Water Resour. Res. 28 (4), 1015-1031.
#'
#' @export
#'
#' @examples
#' set.seed(1234)
#' # Two-dimensional sphere
#' result <- sceua(
#'   fn = function(x) sum(x^2),
#'   lower = c(-5, -5),
#'   upper = c(5, 5),
#'   max_evaluations = 5000,
#'   kstop = 5,
#'   pcento = 1e-8,
#'   complexes = 5
#' )
#' result
sceua <- function(
  fn,
  lower,
  upper,
  initial = NULL,
  max_evaluations = 10000L,
  kstop = 5L,
  pcento = 0.01,
  complexes = 2L,
  points_per_complex = NULL,
  simplex_size = NULL,
  evolution_steps = NULL,
  min_complexes = NULL,
  parameter_epsilon = 1e-3,
  ...
) {
  checkmate::assert_function(fn)
  checkmate::assert_numeric(lower, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(upper, finite = TRUE, any.missing = FALSE)

  if (length(lower) != length(upper)) {
    stop("`lower` and `upper` must have the same length.")
  }
  if (length(lower) == 0L) {
    stop("At least one parameter is required.")
  }
  if (!is.null(initial)) {
    checkmate::assert_numeric(
      initial,
      finite = TRUE,
      any.missing = FALSE,
      len = length(lower)
    )
  }

  n <- length(lower)
  if (is.null(points_per_complex)) {
    points_per_complex <- 2L * n + 1L
  }
  if (is.null(simplex_size)) {
    simplex_size <- n + 1L
  }
  if (is.null(evolution_steps)) {
    evolution_steps <- points_per_complex
  }
  if (is.null(min_complexes)) {
    min_complexes <- complexes
  }

  if (length(list(...)) > 0L) {
    args <- list(...)
    original_fn <- fn
    fn <- function(x) do.call(original_fn, c(list(x), args))
  }

  result <- sceua_minimize(
    fn = fn,
    lower = as.double(lower),
    upper = as.double(upper),
    max_evaluations = as.integer(max_evaluations),
    kstop = as.integer(kstop),
    pcento = as.double(pcento),
    seed = draw_sceua_seed(),
    complexes = as.integer(complexes),
    points_per_complex = as.integer(points_per_complex),
    simplex_size = as.integer(simplex_size),
    evolution_steps = as.integer(evolution_steps),
    min_complexes = as.integer(min_complexes),
    include_initial = !is.null(initial),
    initial_point = if (is.null(initial)) NULL else as.double(initial),
    parameter_epsilon = as.double(parameter_epsilon)
  )

  result$par <- stats::setNames(result$par, names(lower))
  history <- data.frame(
    loop_index = result$history$loop_index,
    evaluations = result$history$evaluations,
    complexes = result$history$complexes,
    best_f = result$history$best_f,
    worst_f = result$history$worst_f,
    geometric_range = result$history$geometric_range,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  history$best_x <- result$history$best_x
  result$history <- history

  structure(result, class = "sceua")
}

draw_sceua_seed <- function() {
  as.integer(stats::runif(1L, min = 1, max = .Machine$integer.max))
}

#' @export
print.sceua <- function(x, ...) {
  cat(
    "<sceua>\n",
    sprintf("best value:    %g\n", x$value),
    sprintf("evaluations:   %d\n", x$counts),
    sprintf("iterations:    %d\n", x$iterations),
    sprintf("termination:   %s\n", x$termination),
    sep = ""
  )
  cat("best parameters:\n")
  print(x$par)
  invisible(x)
}
