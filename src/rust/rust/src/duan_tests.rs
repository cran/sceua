use crate::{duan_test_func::*, minimize, Config};

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

// These optimum checks are direct copies of the documented Matlab test cases by Duan
// See src/duan_test_func.rs for the function definitions
// https://www.mathworks.com/matlabcentral/fileexchange/7671-shuffled-complex-evolution-sce-ua-method

#[test]
fn duan_test_functions_match_documented_optima() {
    assert_close(goldstein_price(&[0.0, -1.0]), 3.0, 1.0e-12);
    assert_close(rosenbrock(&[1.0, 1.0]), 0.0, 1.0e-12);
    assert_close(
        six_hump_camelback(&[0.08983, -0.7126]),
        -1.0316284229280819,
        1.0e-8,
    );
    assert_close(rastrigin_duan(&[0.0, 0.0]), -2.0, 1.0e-12);
    assert_close(griewank_duan(&[0.0, 0.0]), 0.0, 1.0e-12);
    assert_close(griewank_duan(&[0.0; 10]), 0.0, 1.0e-12);
    assert_close(shekel(&[4.0, 4.0, 4.0, 4.0]), -10.536283726219603, 1.0e-12);
    assert_close(
        hartman(&[0.201, 0.150, 0.477, 0.275, 0.311, 0.657]),
        -3.3223349676854577,
        1.0e-12,
    );
}

#[test]
fn sceua_minimizes_goldstein_price_with_duan_bounds() {
    let config = Config {
        max_evaluations: 10_000,
        kstop: 5,
        pcento: 0.01,
        seed: 1969,
        complexes: 5,
        ..Config::default()
    };
    let result = minimize(goldstein_price, &[-2.0, -2.0], &[2.0, 2.0], config).unwrap();

    assert!(result.best_f <= 3.001, "{result:?}");
    assert!((result.best_x[0] - 0.0).abs() <= 0.01, "{result:?}");
    assert!((result.best_x[1] + 1.0).abs() <= 0.01, "{result:?}");
}

#[test]
fn sceua_minimizes_rosenbrock_with_duan_bounds() {
    let config = Config {
        max_evaluations: 20_000,
        kstop: 5,
        pcento: 0.0,
        seed: 1969,
        complexes: 10,
        ..Config::default()
    };
    let result = minimize(rosenbrock, &[-5.0, -5.0], &[5.0, 5.0], config).unwrap();

    assert!(result.best_f <= 1.0e-3, "{result:?}");
    assert!((result.best_x[0] - 1.0).abs() <= 0.05, "{result:?}");
    assert!((result.best_x[1] - 1.0).abs() <= 0.05, "{result:?}");
}
