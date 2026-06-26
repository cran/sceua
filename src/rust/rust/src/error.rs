use std::{error::Error, fmt};

/// Errors returned by the SCE-UA optimizer.
#[derive(Debug, Clone, PartialEq)]
pub enum SceuaError {
    /// No parameters were supplied for optimization.
    EmptyProblem,
    /// Lower and upper parameter-bound vectors have different lengths.
    BoundsLengthMismatch {
        /// Number of lower bounds supplied.
        lower: usize,
        /// Number of upper bounds supplied.
        upper: usize,
    },
    /// Initial parameter set length does not match the number of parameters.
    InitialPointLengthMismatch {
        /// Expected number of parameters.
        expected: usize,
        /// Actual number of values in the initial point.
        actual: usize,
    },
    /// Lower or upper bound on a parameter is invalid.
    InvalidBounds {
        /// Zero-based parameter index.
        index: usize,
        /// Lower bound on the parameter.
        lower: f64,
        /// Upper bound on the parameter.
        upper: f64,
    },
    /// Algorithmic control or convergence-check parameter is invalid.
    InvalidConfig(&'static str),
    /// Objective function returned a non-finite criterion value.
    NonFiniteObjective {
        /// Non-finite value returned by the objective function.
        value: f64,
    },
}

impl fmt::Display for SceuaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyProblem => write!(f, "at least one parameter is required"),
            Self::BoundsLengthMismatch { lower, upper } => write!(
                f,
                "lower and upper bounds must have the same length, got {lower} and {upper}"
            ),
            Self::InitialPointLengthMismatch { expected, actual } => {
                write!(f, "initial point has length {actual}, expected {expected}")
            }
            Self::InvalidBounds {
                index,
                lower,
                upper,
            } => write!(
                f,
                "invalid bounds at index {index}: lower={lower}, upper={upper}"
            ),
            Self::InvalidConfig(message) => write!(f, "invalid SCE-UA configuration: {message}"),
            Self::NonFiniteObjective { value } => {
                write!(f, "objective returned a non-finite value: {value}")
            }
        }
    }
}

impl Error for SceuaError {}
