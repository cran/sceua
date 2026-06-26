use crate::error::SceuaError;

/// SCE-UA algorithmic control and convergence-check parameters.
///
/// The field names mirror the input argument variables in Duan's Fortran
/// implementation where practical.
#[derive(Debug, Clone)]
pub struct Config {
    /// Maximum number of trials allowed before optimization is terminated.
    pub max_evaluations: usize,
    /// Number of shuffling loops in which the criterion value must change by
    /// the given percentage before optimization is terminated.
    pub kstop: usize,
    /// Fractional percentage by which the criterion value must change in
    /// `kstop` shuffling loops.
    ///
    /// For example, `0.01` means one percent.
    pub pcento: f64,
    /// Initial random seed.
    pub seed: i64,
    /// Number of complexes in the initial population.
    pub complexes: usize,
    /// Number of points in each complex.
    ///
    /// Defaults to `2 * n + 1`, where `n` is the number of parameters.
    pub points_per_complex: Option<usize>,
    /// Number of points in a sub-complex.
    ///
    /// Defaults to `n + 1`, where `n` is the number of parameters.
    pub simplex_size: Option<usize>,
    /// Number of evolution steps allowed for each complex before complex
    /// shuffling.
    ///
    /// Defaults to `points_per_complex`.
    pub evolution_steps: Option<usize>,
    /// Minimum number of complexes required if the number of complexes is
    /// allowed to reduce as the optimization proceeds.
    ///
    /// Defaults to `complexes`.
    pub min_complexes: Option<usize>,
    /// Flag on whether to include the initial point in the population.
    pub include_initial: bool,
    /// Initial parameter set.
    pub initial_point: Option<Vec<f64>>,
    /// Parameter convergence threshold for the normalized geometric mean of
    /// parameter ranges.
    pub parameter_epsilon: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            max_evaluations: 10_000,
            kstop: 5,
            pcento: 0.01,
            seed: 1969,
            complexes: 2,
            points_per_complex: None,
            simplex_size: None,
            evolution_steps: None,
            min_complexes: None,
            include_initial: false,
            initial_point: None,
            parameter_epsilon: 1.0e-3,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ResolvedConfig {
    pub max_evaluations: usize,
    pub kstop: usize,
    pub pcento: f64,
    pub seed: i64,
    pub complexes: usize,
    pub points_per_complex: usize,
    pub simplex_size: usize,
    pub evolution_steps: usize,
    pub min_complexes: usize,
    pub parameter_epsilon: f64,
}

impl Config {
    pub(crate) fn resolve(&self, dimension: usize) -> Result<ResolvedConfig, SceuaError> {
        if dimension == 0 {
            return Err(SceuaError::EmptyProblem);
        }
        if self.max_evaluations == 0 {
            return Err(SceuaError::InvalidConfig(
                "max_evaluations must be greater than zero",
            ));
        }
        if self.kstop == 0 || self.kstop > 9 {
            return Err(SceuaError::InvalidConfig(
                "kstop must be in 1..=9 for Duan Fortran compatibility",
            ));
        }
        if !self.pcento.is_finite() || self.pcento < 0.0 {
            return Err(SceuaError::InvalidConfig(
                "pcento must be finite and non-negative",
            ));
        }
        if self.complexes == 0 {
            return Err(SceuaError::InvalidConfig(
                "complexes must be greater than zero",
            ));
        }
        if !self.parameter_epsilon.is_finite() || self.parameter_epsilon <= 0.0 {
            return Err(SceuaError::InvalidConfig(
                "parameter_epsilon must be finite and positive",
            ));
        }

        let points_per_complex = self.points_per_complex.unwrap_or(2 * dimension + 1);
        let simplex_size = self.simplex_size.unwrap_or(dimension + 1);
        let evolution_steps = self.evolution_steps.unwrap_or(points_per_complex);
        let min_complexes = self.min_complexes.unwrap_or(self.complexes);

        if points_per_complex < 2 {
            return Err(SceuaError::InvalidConfig(
                "points_per_complex must be at least 2",
            ));
        }
        if simplex_size < 2 || simplex_size > points_per_complex {
            return Err(SceuaError::InvalidConfig(
                "simplex_size must be in 2..=points_per_complex",
            ));
        }
        if evolution_steps == 0 {
            return Err(SceuaError::InvalidConfig(
                "evolution_steps must be greater than zero",
            ));
        }
        if min_complexes == 0 || min_complexes > self.complexes {
            return Err(SceuaError::InvalidConfig(
                "min_complexes must be in 1..=complexes",
            ));
        }

        Ok(ResolvedConfig {
            max_evaluations: self.max_evaluations,
            kstop: self.kstop,
            pcento: self.pcento,
            seed: self.seed,
            complexes: self.complexes,
            points_per_complex,
            simplex_size,
            evolution_steps,
            min_complexes,
            parameter_epsilon: self.parameter_epsilon,
        })
    }
}
