# sceua

<!-- badges: start -->
[![R-CMD-check](https://github.com/atsyplenkov/sceua/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/atsyplenkov/sceua/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`sceua` provides R bindings to a Rust implementation of the Shuffled Complex
Evolution - University of Arizona (SCE-UA) global optimisation algorithm
(Duan et al., 1992).

SCE-UA combines deterministic simplex search, competitive evolution, and periodic
shuffling of parallel complexes. It is designed for nonlinear, non-convex,
continuous parameter estimation problems and is widely used in hydrological model
calibration.

## Installation

Install the released version from CRAN:

``` r
install.packages("sceua")
```

Or the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("atsyplenkov/sceua/r")
```

Building from source requires a [Rust toolchain](https://rust-lang.org/tools/install/).

## Example

Minimise a simple sphere function:

``` r
library(sceua)
set.seed(1969)
result <- sceua(
  fn = function(x) sum(x^2),
  lower = c(-5, -5),
  upper = c(5, 5),
  max_evaluations = 5000,
  kstop = 5,
  pcento = 1e-8,
  complexes = 5
)

result
#> <sceua>
#> best value:    4.13437e-11
#> evaluations:   625
#> iterations:    15
#> termination:   parameter_convergence
#> best parameters:
#> [1] -5.857289e-06 -2.652528e-06
```

Pass extra arguments to the objective:

``` r
fn <- function(x, target) sum((x - target)^2)
set.seed(1969)
result <- sceua(
  fn = fn,
  lower = c(-5, -5),
  upper = c(5, 5),
  target = c(1, 2),
  max_evaluations = 5000
)

result$par
#> [1] 1.000273 1.999246
```

## Algorithm parameters
The most commonly tuned parameters are:

- `max_evaluations`: maximum number of objective evaluations.
- `kstop`: number of shuffling loops over which the objective must change by
  `pcento` to continue.
- `pcento`: objective convergence threshold (%).
- `complexes`: number of complexes in the initial population.
- `points_per_complex`: points per complex (defaults to `2 * n + 1`).
- `simplex_size`: points in each sub-complex (defaults to `n + 1`).
- `evolution_steps`: evolution steps per complex before shuffling (defaults to
  `points_per_complex`).
- `min_complexes`: minimum number of complexes after reduction (defaults to
  `complexes`).
- `parameter_epsilon`: parameter-space convergence threshold.

See `?sceua` for full details.

## References

Duan, Q., Sorooshian, S., and Gupta, V.K., 1992. Effective and efficient
global optimization for conceptual rainfall-runoff models.
*Water Resources Research* 28 (4), 1015-1031.
