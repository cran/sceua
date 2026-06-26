//! Shuffled Complex Evolution (SCE-UA) method for global optimization.
//!
//! This is a Rust implementation of Duan's (1992) general-purpose global optimization
//! program developed at the Department of Hydrology & Water Resources of the
//! University of Arizona. SCE-UA searches a continuous bounded parameter space
//! by evolving shuffled complexes of points.
//!
//! # Example
//!
//! ```
//! use sceua::{minimize, Config};
//!
//! let result = minimize(
//!     |x| x.iter().map(|value| value * value).sum::<f64>(),
//!     &[-5.0, -5.0],
//!     &[5.0, 5.0],
//!     Config::default(),
//! )?;
//!
//! assert_eq!(result.best_x.len(), 2);
//! # Ok::<(), sceua::SceuaError>(())
//! ```
//!
//! Reference: Duan, Q., Sorooshian, S., and Gupta, V.K. (1992),
//! "Effective and efficient global optimization for conceptual rainfall-runoff
//! models", *Water Resources Research*, 28(4), 1015-1031.

#![warn(missing_docs)]

mod cce;
mod config;
pub mod duan_test_func;
mod error;
mod population;
mod rng;
mod sce;

#[cfg(test)]
mod duan_tests;

pub use config::Config;
pub use error::SceuaError;
pub use sce::{minimize, HistoryEntry, OptimizationResult, TerminationReason};
