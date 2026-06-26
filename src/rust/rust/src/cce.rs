use crate::{error::SceuaError, population::Point, rng::DuanRng};

#[derive(Clone, Copy, Debug)]
pub(crate) struct ComplexLayout {
    offset: usize,
    stride: usize,
}

impl ComplexLayout {
    pub(crate) fn new(offset: usize, stride: usize) -> Self {
        Self { offset, stride }
    }

    fn point_index(self, rank: usize) -> usize {
        rank * self.stride + self.offset
    }
}

#[derive(Debug, Default)]
pub(crate) struct CceScratch {
    worst: Vec<f64>,
    centroid: Vec<f64>,
    step: Vec<f64>,
    trial: Vec<f64>,
}

impl CceScratch {
    fn resize(&mut self, dimension: usize) {
        self.worst.resize(dimension, 0.0);
        self.centroid.resize(dimension, 0.0);
        self.step.resize(dimension, 0.0);
        self.trial.resize(dimension, 0.0);
    }
}

#[derive(Clone, Copy)]
struct SearchSpace<'a> {
    lower: &'a [f64],
    upper: &'a [f64],
    normalized_std: &'a [f64],
}

// CCE subroutine.
// reflection -> contraction -> Gaussian mutation sequence
// See
// https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L431-L546

#[expect(
    clippy::too_many_arguments,
    reason = "CCE mirrors Duan's Fortran/Matlab argument list"
)]
pub(crate) fn evolve_simplex<F>(
    points: &mut [Point],
    layout: ComplexLayout,
    simplex_indices: &[usize],
    lower: &[f64],
    upper: &[f64],
    normalized_std: &[f64],
    scratch: &mut CceScratch,
    rng: &mut DuanRng,
    evaluations: &mut usize,
    max_evaluations: usize,
    objective: &mut F,
) -> Result<Option<usize>, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    if simplex_indices.is_empty() || *evaluations >= max_evaluations {
        return Ok(None);
    }

    let search = SearchSpace {
        lower,
        upper,
        normalized_std,
    };
    let dimension = lower.len();
    scratch.resize(dimension);

    let worst_rank = *simplex_indices
        .last()
        .expect("simplex_indices is not empty");
    let worst_index = layout.point_index(worst_rank);
    scratch.worst.copy_from_slice(&points[worst_index].x);
    let worst_value = points[worst_index].value;

    centroid_without_worst(points, layout, simplex_indices, dimension, scratch);
    for parameter in 0..dimension {
        scratch.step[parameter] = scratch.centroid[parameter] - scratch.worst[parameter];
        scratch.trial[parameter] = scratch.worst[parameter] + 2.0 * scratch.step[parameter];
    }

    if !within_bounds(&scratch.trial, search) {
        gaussian_point(
            points,
            layout,
            simplex_indices[0],
            search,
            rng,
            &mut scratch.trial,
        );
    }

    let mut trial_value = evaluate(objective, &scratch.trial, evaluations)?;
    if trial_value <= worst_value {
        replace_point(&mut points[worst_index], &scratch.trial, trial_value);
        return Ok(Some(worst_rank));
    }
    if *evaluations >= max_evaluations {
        return Ok(None);
    }

    for parameter in 0..dimension {
        scratch.trial[parameter] = scratch.worst[parameter] + 0.5 * scratch.step[parameter];
    }
    trial_value = evaluate(objective, &scratch.trial, evaluations)?;
    if trial_value <= worst_value {
        replace_point(&mut points[worst_index], &scratch.trial, trial_value);
        return Ok(Some(worst_rank));
    }
    if *evaluations >= max_evaluations {
        return Ok(None);
    }

    gaussian_point(
        points,
        layout,
        simplex_indices[0],
        search,
        rng,
        &mut scratch.trial,
    );
    trial_value = evaluate(objective, &scratch.trial, evaluations)?;
    replace_point(&mut points[worst_index], &scratch.trial, trial_value);
    Ok(Some(worst_rank))
}

fn centroid_without_worst(
    points: &[Point],
    layout: ComplexLayout,
    simplex_indices: &[usize],
    dimension: usize,
    scratch: &mut CceScratch,
) {
    let divisor = (simplex_indices.len() - 1) as f64;
    scratch.centroid.fill(0.0);
    for &rank in &simplex_indices[..simplex_indices.len() - 1] {
        let point = &points[layout.point_index(rank)];
        for (sum, value) in scratch.centroid.iter_mut().zip(&point.x) {
            *sum += *value;
        }
    }
    let inv = 1.0 / divisor;
    for value in scratch.centroid.iter_mut().take(dimension) {
        *value *= inv;
    }
}

fn gaussian_point(
    points: &[Point],
    layout: ComplexLayout,
    best_rank: usize,
    search: SearchSpace<'_>,
    rng: &mut DuanRng,
    trial: &mut [f64],
) {
    let best = &points[layout.point_index(best_rank)].x;
    for (parameter, value) in trial.iter_mut().enumerate() {
        let bound = search.upper[parameter] - search.lower[parameter];
        loop {
            let candidate =
                best[parameter] + search.normalized_std[parameter] * rng.gaussian() * bound;
            if candidate >= search.lower[parameter] && candidate <= search.upper[parameter] {
                *value = candidate;
                break;
            }
        }
    }
}

fn replace_point(point: &mut Point, trial: &[f64], value: f64) {
    point.x.clear();
    point.x.extend_from_slice(trial);
    point.value = value;
}

fn within_bounds(point: &[f64], search: SearchSpace<'_>) -> bool {
    point
        .iter()
        .zip(search.lower.iter().zip(search.upper))
        .all(|(&value, (&lo, &hi))| value >= lo && value <= hi)
}

fn evaluate<F>(objective: &mut F, point: &[f64], evaluations: &mut usize) -> Result<f64, SceuaError>
where
    F: FnMut(&[f64]) -> f64,
{
    let value = objective(point);
    *evaluations += 1;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(SceuaError::NonFiniteObjective { value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: f64, value: f64) -> Point {
        Point { x: vec![x], value }
    }

    fn layout() -> ComplexLayout {
        ComplexLayout::new(0, 1)
    }

    // Mirrors Fortran CCE paths: reflection, contraction, mutation, and maxn stop.
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L431-L546

    #[test]
    fn cce_accepts_reflection_when_it_improves_worst_point() {
        let mut simplex = vec![point(0.0, 0.0), point(1.0, 1.0), point(2.0, 4.0)];
        let mut rng = DuanRng::new(1969);
        let mut evaluations = 0;
        let mut objective = |x: &[f64]| x[0] * x[0];

        let mut scratch = CceScratch::default();
        let changed = evolve_simplex(
            &mut simplex,
            layout(),
            &[0, 1, 2],
            &[-10.0],
            &[10.0],
            &[0.1],
            &mut scratch,
            &mut rng,
            &mut evaluations,
            10,
            &mut objective,
        )
        .unwrap();

        assert_eq!(changed, Some(2));

        assert_eq!(evaluations, 1);
        assert!((simplex[2].x[0] + 1.0).abs() < 1.0e-12);
        assert!((simplex[2].value - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn cce_tries_contraction_after_failed_reflection() {
        let mut simplex = vec![point(0.0, 0.0), point(1.0, 0.5), point(2.0, 1.0)];
        let mut rng = DuanRng::new(1969);
        let mut evaluations = 0;
        let mut objective = |x: &[f64]| (x[0] - 1.0) * (x[0] - 1.0);

        let mut scratch = CceScratch::default();
        let changed = evolve_simplex(
            &mut simplex,
            layout(),
            &[0, 1, 2],
            &[-10.0],
            &[10.0],
            &[0.1],
            &mut scratch,
            &mut rng,
            &mut evaluations,
            10,
            &mut objective,
        )
        .unwrap();

        assert_eq!(changed, Some(2));

        assert_eq!(evaluations, 2);
        assert!((simplex[2].x[0] - 1.25).abs() < 1.0e-12);
        assert!((simplex[2].value - 0.0625).abs() < 1.0e-12);
    }

    #[test]
    fn cce_uses_gaussian_mutation_when_reflection_is_out_of_bounds() {
        let mut simplex = vec![point(0.9, 0.0), point(0.1, 1.0)];
        let mut rng = DuanRng::new(1969);
        let mut evaluations = 0;
        let mut objective = |x: &[f64]| (x[0] - 0.9).abs();

        let mut scratch = CceScratch::default();
        let changed = evolve_simplex(
            &mut simplex,
            layout(),
            &[0, 1],
            &[0.0],
            &[1.0],
            &[0.05],
            &mut scratch,
            &mut rng,
            &mut evaluations,
            10,
            &mut objective,
        )
        .unwrap();

        assert_eq!(changed, Some(1));

        assert_eq!(evaluations, 1);
        assert!(simplex[1].x[0] >= 0.0 && simplex[1].x[0] <= 1.0);
        assert_ne!(simplex[1].x[0], 1.7);
    }

    #[test]
    fn cce_does_not_replace_worst_when_failed_reflection_hits_max_evaluations() {
        let mut simplex = vec![point(0.0, 0.0), point(1.0, 0.5), point(2.0, 1.0)];
        let mut rng = DuanRng::new(1969);
        let mut evaluations = 0;
        let mut objective = |_x: &[f64]| 2.0;

        let mut scratch = CceScratch::default();
        let changed = evolve_simplex(
            &mut simplex,
            layout(),
            &[0, 1, 2],
            &[-10.0],
            &[10.0],
            &[0.1],
            &mut scratch,
            &mut rng,
            &mut evaluations,
            1,
            &mut objective,
        )
        .unwrap();

        assert_eq!(changed, None);

        assert_eq!(evaluations, 1);
        assert_eq!(simplex[2], point(2.0, 1.0));
    }
}
