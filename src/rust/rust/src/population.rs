use crate::rng::DuanRng;

const DELTA: f64 = 1.0e-20;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Point {
    pub(crate) x: Vec<f64>,
    pub(crate) value: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ParameterStats {
    pub(crate) normalized_std: Vec<f64>,
    pub(crate) geometric_range: f64,
}

pub(crate) fn sort_points(points: &mut [Point]) {
    points.sort_by(|left, right| left.value.total_cmp(&right.value));
}

pub(crate) fn random_point(lower: &[f64], upper: &[f64], rng: &mut DuanRng) -> Vec<f64> {
    lower
        .iter()
        .zip(upper)
        .map(|(&lo, &hi)| lo + (hi - lo) * rng.uniform())
        .collect()
}

// PARSTT computes normalised population std and geometric range. See
// https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L552-L589

pub(crate) fn parameter_stats(points: &[Point], lower: &[f64], upper: &[f64]) -> ParameterStats {
    let dimension = lower.len();
    let count = points.len() as f64;
    let mut normalized_std = vec![0.0; dimension];
    let mut log_range_sum = 0.0;

    for parameter in 0..dimension {
        let bound = upper[parameter] - lower[parameter];
        let mut min_value = f64::INFINITY;
        let mut max_value = f64::NEG_INFINITY;
        let mut sum = 0.0;
        let mut sum_squares = 0.0;

        for point in points {
            debug_assert_eq!(point.x.len(), dimension);
            let value = point.x[parameter];
            min_value = min_value.min(value);
            max_value = max_value.max(value);
            sum += value;
            sum_squares += value * value;
        }

        let mean = sum / count;
        let mut variance = sum_squares / count - mean * mean;
        if variance <= DELTA {
            variance = DELTA;
        }
        normalized_std[parameter] = variance.sqrt() / bound;
        log_range_sum += (DELTA + (max_value - min_value) / bound).ln();
    }

    ParameterStats {
        normalized_std,
        geometric_range: (log_range_sum / dimension as f64).exp(),
    }
}

#[cfg(test)]
pub(crate) fn normalized_distances(
    points: &[Point],
    initial: &[f64],
    lower: &[f64],
    upper: &[f64],
) -> Vec<f64> {
    points
        .iter()
        .map(|point| {
            point
                .x
                .iter()
                .zip(initial)
                .zip(lower.iter().zip(upper))
                .map(|((&x, &xi), (&lo, &hi))| (x - xi).abs() / (hi - lo))
                .sum::<f64>()
                / initial.len() as f64
        })
        .collect()
}

// linear-probability sub-complex sampling before SORT1.
// See https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L240-L267

#[derive(Debug, Default)]
pub(crate) struct SimplexSampler {
    indices: Vec<usize>,
    marks: Vec<u32>,
    generation: u32,
}

impl SimplexSampler {
    pub(crate) fn sample(
        &mut self,
        points_per_complex: usize,
        simplex_size: usize,
        rng: &mut DuanRng,
    ) -> &[usize] {
        self.indices.clear();

        if simplex_size == points_per_complex {
            self.indices.extend(0..simplex_size);
            return &self.indices;
        }

        if self.marks.len() < points_per_complex {
            self.marks.resize(points_per_complex, 0);
        }
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.marks.fill(0);
            self.generation = 1;
        }

        while self.indices.len() < simplex_size {
            let candidate = sample_rank(points_per_complex, rng);
            if self.marks[candidate] != self.generation {
                self.marks[candidate] = self.generation;
                self.indices.push(candidate);
            }
        }
        self.indices.sort_unstable();
        &self.indices
    }
}

#[cfg(test)]
pub(crate) fn sample_simplex_indices(
    points_per_complex: usize,
    simplex_size: usize,
    rng: &mut DuanRng,
) -> Vec<usize> {
    let mut sampler = SimplexSampler::default();
    sampler
        .sample(points_per_complex, simplex_size, rng)
        .to_vec()
}

// COMP drops the lowest-ranked complex during reduction.
// See https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L618-L644

pub(crate) fn compress_complexes(
    population: Vec<Point>,
    old_complexes: usize,
    new_complexes: usize,
    points_per_complex: usize,
) -> Vec<Point> {
    let mut compressed = Vec::with_capacity(new_complexes * points_per_complex);
    let mut population = population.into_iter();
    for _ in 0..points_per_complex {
        for complex_index in 0..old_complexes {
            let point = population.next().expect("population has expected size");
            if complex_index < new_complexes {
                compressed.push(point);
            }
        }
    }
    compressed
}
fn sample_rank(points_per_complex: usize, rng: &mut DuanRng) -> usize {
    let npg = points_per_complex as f64;
    let random = rng.uniform();
    let npg_half = npg + 0.5;
    let one_based =
        1.0 + (npg_half - (npg_half * npg_half - npg * (npg + 1.0) * random).sqrt()).trunc();
    one_based as usize - 1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: &[f64], value: f64) -> Point {
        Point {
            x: x.to_vec(),
            value,
        }
    }

    // Mirrors Fortran PARSTT normalised std plus normalised geometric range
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L552-L589

    #[test]
    fn parstt_matches_duan_population_statistics() {
        let points = vec![
            point(&[0.0, 2.0], 0.0),
            point(&[2.0, 4.0], 0.0),
            point(&[4.0, 8.0], 0.0),
        ];
        let stats = parameter_stats(&points, &[0.0, 0.0], &[4.0, 8.0]);

        let expected_std = (8.0_f64 / 3.0).sqrt() / 4.0;
        assert!((stats.normalized_std[0] - expected_std).abs() < 1.0e-12);
        assert!((stats.normalized_std[1] - 0.31180478223116176).abs() < 1.0e-12);
        assert!((stats.geometric_range - 0.75_f64.sqrt()).abs() < 1.0e-12);
    }

    // Mirrors Fortran NORMDIST; used for reported population distance.
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L594-L613

    #[test]
    fn normdist_matches_duan_formula() {
        let points = vec![point(&[2.0, 6.0], 0.0), point(&[4.0, 0.0], 0.0)];
        let distances = normalized_distances(&points, &[0.0, 4.0], &[0.0, 0.0], &[4.0, 8.0]);
        assert!((distances[0] - 0.375).abs() < 1.0e-12);
        assert!((distances[1] - 0.75).abs() < 1.0e-12);
    }

    // Mirrors Fortran sub-complex selection formula and SORT1 ordering.
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L240-L267

    #[test]
    fn simplex_indices_follow_duan_linear_probability() {
        let mut rng = DuanRng::new(1969);
        let indices = sample_simplex_indices(5, 3, &mut rng);
        assert_eq!(indices, vec![0, 2, 3]);
    }

    // Mirrors Fortran COMP complex-number reduction.
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L618-L644

    #[test]
    fn comp_drops_lowest_ranked_complex() {
        let population: Vec<_> = (0..9).map(|i| point(&[i as f64], i as f64)).collect();
        let compressed = compress_complexes(population, 3, 2, 3);
        let values: Vec<_> = compressed
            .iter()
            .map(|point| point.value as usize)
            .collect();
        assert_eq!(values, vec![0, 1, 3, 4, 6, 7]);
    }
}
