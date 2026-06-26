// NB: RAN1/GASDEV in FUSE SCE, kept in single precision to
// match the original REAL arithmetic. See fortran source:
// https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L788-L849

#[derive(Clone, Debug)]
pub(crate) struct DuanRng {
    idum: i64,
    initialized: bool,
    ix1: i64,
    ix2: i64,
    ix3: i64,
    table: [f32; 97],
    has_gaussian: bool,
    gaussian: f32,
}

impl DuanRng {
    pub(crate) fn new(seed: i64) -> Self {
        let idum = seed.checked_abs().map_or(-i64::MAX, |seed| -seed);
        Self {
            idum,
            initialized: false,
            ix1: 0,
            ix2: 0,
            ix3: 0,
            table: [0.0; 97],
            has_gaussian: false,
            gaussian: 0.0,
        }
    }

    pub(crate) fn uniform(&mut self) -> f64 {
        f64::from(self.raw_uniform())
    }

    fn raw_uniform(&mut self) -> f32 {
        const M1: i64 = 259_200;
        const IA1: i64 = 7_141;
        const IC1: i64 = 54_773;
        const RM1: f32 = 3.858_024_7e-6;
        const M2: i64 = 134_456;
        const IA2: i64 = 8_121;
        const IC2: i64 = 28_411;
        const RM2: f32 = 7.437_377_3e-6;
        const M3: i64 = 243_000;
        const IA3: i64 = 4_561;
        const IC3: i64 = 51_349;

        if self.idum < 0 || !self.initialized {
            self.initialized = true;
            self.ix1 = (IC1 - self.idum).rem_euclid(M1);
            self.ix1 = ((IA1 * self.ix1) + IC1).rem_euclid(M1);
            self.ix2 = self.ix1.rem_euclid(M2);
            self.ix1 = ((IA1 * self.ix1) + IC1).rem_euclid(M1);
            self.ix3 = self.ix1.rem_euclid(M3);
            for value in &mut self.table {
                self.ix1 = ((IA1 * self.ix1) + IC1).rem_euclid(M1);
                self.ix2 = ((IA2 * self.ix2) + IC2).rem_euclid(M2);
                *value = ran1_value(self.ix1, self.ix2, RM1, RM2);
            }
            self.idum = 1;
            self.has_gaussian = false;
        }

        self.ix1 = ((IA1 * self.ix1) + IC1).rem_euclid(M1);
        self.ix2 = ((IA2 * self.ix2) + IC2).rem_euclid(M2);
        self.ix3 = ((IA3 * self.ix3) + IC3).rem_euclid(M3);
        let j = ((97 * self.ix3) / M3) as usize;
        let random = self.table[j];
        self.table[j] = ran1_value(self.ix1, self.ix2, RM1, RM2);
        random
    }

    pub(crate) fn gaussian(&mut self) -> f64 {
        if self.has_gaussian {
            self.has_gaussian = false;
            return f64::from(self.gaussian);
        }

        loop {
            let v1 = 2.0_f32 * self.raw_uniform() - 1.0;
            let v2 = 2.0_f32 * self.raw_uniform() - 1.0;
            let radius = v1 * v1 + v2 * v2;
            if radius > 0.0 && radius < 1.0 {
                let factor = (-(2.0_f32 * radius.ln() / radius)).sqrt();
                self.gaussian = v1 * factor;
                self.has_gaussian = true;
                return f64::from(v2 * factor);
            }
        }
    }
}

fn ran1_value(ix1: i64, ix2: i64, rm1: f32, rm2: f32) -> f32 {
    (((ix1 as f64) + (ix2 as f64) * f64::from(rm2)) * f64::from(rm1)) as f32
}

#[cfg(test)]
mod tests {
    use super::DuanRng;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-6,
            "actual={actual}, expected={expected}"
        );
    }

    // Golden values generated from Duan/FUSE Fortran RAN1 with seed 1969:
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L788-L824

    #[test]
    fn ran1_matches_duan_fortran_sequence() {
        let expected = [
            0.051818177104,
            0.834282934666,
            0.799921333790,
            0.820423126221,
            0.492580622435,
            0.482416421175,
            0.574328243732,
            0.130261763930,
            0.688977003098,
            0.465039193630,
        ];
        let mut rng = DuanRng::new(1969);
        for expected in expected {
            assert_close(rng.uniform(), expected);
        }
    }

    // Golden values generated from Duan/FUSE Fortran GASDEV with static SAVE state:
    // https://github.com/naddor/fuse/blob/e5fe0fbed82125eec4711854e1c5492da254df41/build/FUSE_SRC/FUSE_SCE/sce.f#L829-L849

    #[test]
    fn gasdev_matches_duan_fortran_static_sequence() {
        let expected = [
            0.527196645737,
            0.493464767933,
            -3.329959630966,
            -1.405074000359,
            -1.041249513626,
            0.209321752191,
        ];
        let mut rng = DuanRng::new(1969);
        for expected in expected {
            assert_close(rng.gaussian(), expected);
        }
    }
}
