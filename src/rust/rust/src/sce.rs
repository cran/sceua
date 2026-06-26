use crate::{
    cce::{evolve_simplex, CceScratch, ComplexLayout},
    config::{Config, ResolvedConfig},
    error::SceuaError,
    population::{
        compress_complexes, parameter_stats, random_point, sort_points, ParameterStats, Point,
        SimplexSampler,
    },
    rng::DuanRng,
};

/// Final parameter estimate, criterion value, and search metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct OptimizationResult {
    /// Best point found by the optimization search.
    pub best_x: Vec<f64>,
    /// Function value of `best_x`.
    pub best_f: f64,
    /// Number of function evaluations used.
    pub evaluations: usize,
    /// Number of completed shuffling loops.
    pub loops: usize,
    /// Reason the search terminated.
    pub termination: TerminationReason,
    /// Results recorded for the initial population and each shuffling loop.
    pub history: Vec<HistoryEntry>,
}

/// Population summary recorded after shuffling loops.
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryEntry {
    /// Shuffling-loop index. The initial population is recorded as loop `0`.
    pub loop_index: usize,
    /// Number of function evaluations used at this loop.
    pub evaluations: usize,
    /// Number of complexes in the current population.
    pub complexes: usize,
    /// Best criterion value in the current population.
    pub best_f: f64,
    /// Worst criterion value in the current population.
    pub worst_f: f64,
    /// Normalized geometric mean of parameter ranges.
    pub geometric_range: f64,
    /// Best point in the current population.
    pub best_x: Vec<f64>,
}

/// Search termination condition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TerminationReason {
    /// Search terminated because the limit on the maximum number of trials was
    /// exceeded.
    MaxEvaluations,
    /// Optimization terminated because the criterion value has not changed by
    /// the configured percentage in the configured number of shuffling loops.
    ObjectiveConvergence,
    /// Optimization terminated because the population has converged into a
    /// sufficiently small percentage of the feasible space.
    ParameterConvergence,
}

/// Run the Shuffled Complex Evolution method for global optimization.
///
/// The objective function receives a parameter set and must return the
/// criterion value to minimize. `lower` and `upper` are the lower and upper
/// bounds on the parameters.
///
/// This follows Duan's SCE-UA main routine: generate an initial population,
/// arrange points in order of increasing function value, evolve complexes,
/// shuffle, and stop when one of the convergence checks is satisfied.
///
/// Source routine: <https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L152-L399>
pub fn minimize<F>(
    mut objective: F,
    lower: &[f64],
    upper: &[f64],
    config: Config,
) -> Result<OptimizationResult, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    validate_problem(lower, upper, &config)?;
    let resolved = config.resolve(lower.len())?;
    let mut rng = DuanRng::new(resolved.seed);
    let (population, evaluations) =
        initialize_population_serial(&mut objective, lower, upper, &config, resolved, &mut rng)?;
    continue_minimize(
        objective,
        lower,
        upper,
        resolved,
        population,
        evaluations,
        rng,
    )
}

fn initialize_population_serial<F>(
    objective: &mut F,
    lower: &[f64],
    upper: &[f64],
    config: &Config,
    resolved: ResolvedConfig,
    rng: &mut DuanRng,
) -> Result<(Vec<Point>, usize), SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    let target_population = resolved.complexes * resolved.points_per_complex;
    let mut evaluations = 0usize;
    let mut population = Vec::with_capacity(target_population);

    if config.include_initial {
        let initial = config
            .initial_point
            .as_deref()
            .ok_or(SceuaError::InvalidConfig(
                "include_initial requires initial_point",
            ))?;
        population.push(evaluate_point(objective, initial, &mut evaluations)?);
    }

    while population.len() < target_population && evaluations < resolved.max_evaluations {
        let point = random_point(lower, upper, rng);
        population.push(evaluate_owned_point(objective, point, &mut evaluations)?);
    }

    Ok((population, evaluations))
}

fn continue_minimize<F>(
    mut objective: F,
    lower: &[f64],
    upper: &[f64],
    resolved: ResolvedConfig,
    mut population: Vec<Point>,
    evaluations: usize,
    mut rng: DuanRng,
) -> Result<OptimizationResult, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    sort_points(&mut population);
    let mut current_complexes = resolved.complexes;
    let mut evaluations = evaluations;
    let mut history = Vec::new();
    let mut current_stats = parameter_stats(&population, lower, upper);
    push_history(
        &mut history,
        0,
        evaluations,
        current_complexes,
        &population,
        &current_stats,
    );

    if evaluations >= resolved.max_evaluations {
        return Ok(result(
            population,
            evaluations,
            0,
            TerminationReason::MaxEvaluations,
            history,
        ));
    }

    if current_stats.geometric_range <= resolved.parameter_epsilon {
        return Ok(result(
            population,
            evaluations,
            0,
            TerminationReason::ParameterConvergence,
            history,
        ));
    }

    let mut best_by_loop = Vec::new();
    let mut sampler = SimplexSampler::default();
    let mut cce_scratch = CceScratch::default();
    let mut loops = 0usize;

    loop {
        loops += 1;

        for complex_index in 0..current_complexes {
            let layout = ComplexLayout::new(complex_index, current_complexes);
            for _ in 0..resolved.evolution_steps {
                if evaluations >= resolved.max_evaluations {
                    break;
                }

                let simplex_indices =
                    sampler.sample(resolved.points_per_complex, resolved.simplex_size, &mut rng);

                if let Some(changed_rank) = evolve_simplex(
                    &mut population,
                    layout,
                    simplex_indices,
                    lower,
                    upper,
                    &current_stats.normalized_std,
                    &mut cce_scratch,
                    &mut rng,
                    &mut evaluations,
                    resolved.max_evaluations,
                    &mut objective,
                )? {
                    reposition_complex_point(
                        &mut population,
                        complex_index,
                        current_complexes,
                        resolved.points_per_complex,
                        changed_rank,
                    );
                }
            }

            if evaluations >= resolved.max_evaluations {
                break;
            }
        }

        sort_points(&mut population);
        let next_stats = parameter_stats(&population, lower, upper);
        push_history(
            &mut history,
            loops,
            evaluations,
            current_complexes,
            &population,
            &next_stats,
        );
        best_by_loop.push(population[0].value);

        if evaluations >= resolved.max_evaluations {
            return Ok(result(
                population,
                evaluations,
                loops,
                TerminationReason::MaxEvaluations,
                history,
            ));
        }

        if best_by_loop.len() > resolved.kstop {
            let current = *best_by_loop.last().expect("best_by_loop is not empty");
            let old = best_by_loop[0];
            let denominator = (old + current).abs() / 2.0;
            let timeout = if denominator == 0.0 {
                if old == current {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                (old - current).abs() / denominator
            };
            if timeout < resolved.pcento {
                return Ok(result(
                    population,
                    evaluations,
                    loops,
                    TerminationReason::ObjectiveConvergence,
                    history,
                ));
            }
            best_by_loop.remove(0);
        }

        if next_stats.geometric_range <= resolved.parameter_epsilon {
            return Ok(result(
                population,
                evaluations,
                loops,
                TerminationReason::ParameterConvergence,
                history,
            ));
        }

        current_stats = next_stats;
        if current_complexes > resolved.min_complexes {
            let reduced = current_complexes - 1;
            population = compress_complexes(
                population,
                current_complexes,
                reduced,
                resolved.points_per_complex,
            );
            current_complexes = reduced;
        }
    }
}

fn validate_problem(lower: &[f64], upper: &[f64], config: &Config) -> Result<(), SceuaError> {
    if lower.len() != upper.len() {
        return Err(SceuaError::BoundsLengthMismatch {
            lower: lower.len(),
            upper: upper.len(),
        });
    }
    if lower.is_empty() {
        return Err(SceuaError::EmptyProblem);
    }
    for (index, (&lo, &hi)) in lower.iter().zip(upper).enumerate() {
        if !lo.is_finite() || !hi.is_finite() || hi <= lo {
            return Err(SceuaError::InvalidBounds {
                index,
                lower: lo,
                upper: hi,
            });
        }
    }
    if let Some(initial) = &config.initial_point {
        if initial.len() != lower.len() {
            return Err(SceuaError::InitialPointLengthMismatch {
                expected: lower.len(),
                actual: initial.len(),
            });
        }
        for (index, ((&value, &lo), &hi)) in initial.iter().zip(lower).zip(upper).enumerate() {
            if !value.is_finite() || value < lo || value > hi {
                return Err(SceuaError::InvalidBounds {
                    index,
                    lower: lo,
                    upper: hi,
                });
            }
        }
    }
    Ok(())
}

fn evaluate_point<F>(
    objective: &mut F,
    point: &[f64],
    evaluations: &mut usize,
) -> Result<Point, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    let value = objective(point);
    *evaluations += 1;
    if value.is_finite() {
        Ok(Point {
            x: point.to_vec(),
            value,
        })
    } else {
        Err(SceuaError::NonFiniteObjective { value })
    }
}

fn evaluate_owned_point<F>(
    objective: &mut F,
    point: Vec<f64>,
    evaluations: &mut usize,
) -> Result<Point, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    let value = objective(&point);
    *evaluations += 1;
    if value.is_finite() {
        Ok(Point { x: point, value })
    } else {
        Err(SceuaError::NonFiniteObjective { value })
    }
}

fn reposition_complex_point(
    population: &mut [Point],
    complex_index: usize,
    complexes: usize,
    points_per_complex: usize,
    mut rank: usize,
) {
    while rank > 0 {
        let current = rank * complexes + complex_index;
        let previous = (rank - 1) * complexes + complex_index;
        if population[current]
            .value
            .total_cmp(&population[previous].value)
            .is_lt()
        {
            population.swap(current, previous);
            rank -= 1;
        } else {
            break;
        }
    }

    while rank + 1 < points_per_complex {
        let current = rank * complexes + complex_index;
        let next = (rank + 1) * complexes + complex_index;
        if population[next]
            .value
            .total_cmp(&population[current].value)
            .is_lt()
        {
            population.swap(current, next);
            rank += 1;
        } else {
            break;
        }
    }
}

fn push_history(
    history: &mut Vec<HistoryEntry>,
    loop_index: usize,
    evaluations: usize,
    complexes: usize,
    population: &[Point],
    stats: &ParameterStats,
) {
    let best = &population[0];
    let worst = population.last().expect("population is not empty");
    history.push(HistoryEntry {
        loop_index,
        evaluations,
        complexes,
        best_f: best.value,
        worst_f: worst.value,
        geometric_range: stats.geometric_range,
        best_x: best.x.clone(),
    });
}

fn result(
    population: Vec<Point>,
    evaluations: usize,
    loops: usize,
    termination: TerminationReason,
    history: Vec<HistoryEntry>,
) -> OptimizationResult {
    let best = &population[0];
    OptimizationResult {
        best_x: best.x.clone(),
        best_f: best.value,
        evaluations,
        loops,
        termination,
        history,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: f64, value: f64) -> Point {
        Point { x: vec![x], value }
    }

    #[test]
    fn reposition_complex_point_matches_stable_sort_order() {
        let mut population = vec![
            point(0.0, 1.0),
            point(1.0, 2.0),
            point(2.0, 2.0),
            point(3.0, 3.0),
        ];
        population[1].value = 3.0;
        let mut expected = population.clone();
        sort_points(&mut expected);

        reposition_complex_point(&mut population, 0, 1, 4, 1);

        assert_eq!(population, expected);
    }

    #[test]
    fn reposition_complex_point_uses_strided_complex_layout() {
        let mut population = vec![
            point(0.0, 0.0),
            point(10.0, 1.0),
            point(1.0, 2.0),
            point(11.0, 5.0),
            point(2.0, 4.0),
            point(12.0, 6.0),
        ];
        population[3].value = 7.0;

        reposition_complex_point(&mut population, 1, 2, 3, 1);

        assert_eq!(population[1], point(10.0, 1.0));
        assert_eq!(population[3], point(12.0, 6.0));
        assert_eq!(population[5], point(11.0, 7.0));
    }

    #[test]
    fn minimize_rejects_mismatched_bounds() {
        let err = minimize(|x| x[0], &[0.0], &[1.0, 2.0], Config::default()).unwrap_err();
        assert_eq!(err, SceuaError::BoundsLengthMismatch { lower: 1, upper: 2 });
    }

    // Full optimiser test identical to the Fortran version.
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L152-L399

    #[test]
    fn minimize_converges_on_two_dimensional_sphere() {
        let config = Config {
            max_evaluations: 5_000,
            kstop: 5,
            pcento: 1.0e-8,
            seed: 1969,
            complexes: 5,
            ..Config::default()
        };
        let result = minimize(
            |x| x.iter().map(|value| value * value).sum::<f64>(),
            &[-5.0, -5.0],
            &[5.0, 5.0],
            config,
        )
        .unwrap();

        assert!(result.best_f < 1.0e-6, "{result:?}");
        assert!(matches!(
            result.termination,
            TerminationReason::ObjectiveConvergence
                | TerminationReason::ParameterConvergence
                | TerminationReason::MaxEvaluations
        ));
    }
}
