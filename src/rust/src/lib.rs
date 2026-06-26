use extendr_api::error::{Error, Result};
use extendr_api::prelude::*;
use sceua::{minimize, Config, OptimizationResult, TerminationReason};
use std::cell::RefCell;
use std::rc::Rc;

#[allow(clippy::too_many_arguments)]
#[extendr]
fn sceua_minimize(
    fun: Robj,
    lower: Vec<f64>,
    upper: Vec<f64>,
    max_evaluations: i32,
    kstop: i32,
    pcento: f64,
    seed: i32,
    complexes: i32,
    points_per_complex: Option<i32>,
    simplex_size: Option<i32>,
    evolution_steps: Option<i32>,
    min_complexes: Option<i32>,
    include_initial: bool,
    initial_point: Option<Vec<f64>>,
    parameter_epsilon: f64,
) -> Result<Robj> {
    let rfun = fun
        .as_function()
        .ok_or_else(|| Error::Other("`fn` must be an R function.".into()))?;

    let call_error: Rc<RefCell<Option<Error>>> = Rc::new(RefCell::new(None));
    let call_error_clone = call_error.clone();
    let mut objective = move |x: &[f64]| -> f64 {
        let x_robj = Robj::from(x.to_vec());
        match rfun.call(pairlist!(x = x_robj)) {
            Ok(value) => value
                .as_real()
                .unwrap_or_else(|| value.as_integer().map_or(f64::NAN, |v| v as f64)),
            Err(error) => {
                if call_error_clone.borrow().is_none() {
                    *call_error_clone.borrow_mut() = Some(error);
                }
                f64::NAN
            }
        }
    };

    let config = Config {
        max_evaluations: max_evaluations as usize,
        kstop: kstop as usize,
        pcento,
        seed: seed as i64,
        complexes: complexes as usize,
        points_per_complex: points_per_complex.map(|v| v as usize),
        simplex_size: simplex_size.map(|v| v as usize),
        evolution_steps: evolution_steps.map(|v| v as usize),
        min_complexes: min_complexes.map(|v| v as usize),
        include_initial,
        initial_point,
        parameter_epsilon,
    };

    let result = minimize(&mut objective, &lower, &upper, config)
        .map_err(|error| Error::Other(error.to_string()))?;

    if let Some(error) = call_error.borrow_mut().take() {
        return Err(error);
    }

    Ok(result_to_robj(result))
}

fn result_to_robj(result: OptimizationResult) -> Robj {
    let loop_index: Vec<i32> = result
        .history
        .iter()
        .map(|entry| entry.loop_index as i32)
        .collect();
    let evaluations: Vec<i32> = result
        .history
        .iter()
        .map(|entry| entry.evaluations as i32)
        .collect();
    let complexes: Vec<i32> = result
        .history
        .iter()
        .map(|entry| entry.complexes as i32)
        .collect();
    let best_f: Vec<f64> = result.history.iter().map(|entry| entry.best_f).collect();
    let worst_f: Vec<f64> = result.history.iter().map(|entry| entry.worst_f).collect();
    let geometric_range: Vec<f64> = result
        .history
        .iter()
        .map(|entry| entry.geometric_range)
        .collect();
    let best_x: Vec<Robj> = result
        .history
        .iter()
        .map(|entry| Robj::from(entry.best_x.clone()))
        .collect();

    let termination = match result.termination {
        TerminationReason::MaxEvaluations => "max_evaluations",
        TerminationReason::ObjectiveConvergence => "objective_convergence",
        TerminationReason::ParameterConvergence => "parameter_convergence",
    };

    list!(
        par = result.best_x,
        value = result.best_f,
        counts = result.evaluations as i32,
        iterations = result.loops as i32,
        termination = termination,
        history = list!(
            loop_index = loop_index,
            evaluations = evaluations,
            complexes = complexes,
            best_f = best_f,
            worst_f = worst_f,
            geometric_range = geometric_range,
            best_x = List::from_values(best_x)
        )
    )
    .into()
}

extendr_module! {
    mod sceua;
    fn sceua_minimize;
}
