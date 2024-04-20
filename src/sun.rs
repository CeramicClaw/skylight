use crate::*;

pub const DEBUG: bool = true;

pub struct Sun {
    pub date: DateTime,
}

impl std::fmt::Display for Sun {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}",
        self.date,
        )
    }
}

impl Sun {
    pub fn new_day() -> Sun {
        Sun {
            date: new_day(2000, Month::JANUARY, 1),
        }
    }

    /// Get a characteristic sun at a given DateTime and observer lat/lon
    pub fn set(&mut self, date: DateTime, obs_lat_deg: f64, obs_lon_deg: f64, obs_elev_m: f64, p_mbar: f64, t_deg_c: f64) {
        self.date = date;
        // Julian dates
        let jd = julian_day(&self.date);
        let jde = julian_ephemeris_day(jd, self.date.delta_t);
        let jc = julian_century(jd);
        let jce = julian_ephemeris_century(jde);
        let jme = jce / 10.0;

        // Helicentric longitude (degrees)
        let l_deg = l_deg(jme); 

        // Heliocentric Latitude (degrees)
        let b_deg = b_deg(jme);
        
        // Earth Radius Vector (AU)
        let r_au = r_au(jme);

        // Geocentric Longitude (degrees)
        let theta_deg = (l_deg + 180.0) % 360.0;

        // Geocentric Latitude (degrees)
        let beta_deg = b_deg * -1.0;

        // Mean elogation of the moon from the sun (degrees)
        let x0 = 297.85036 + 445267.111480 * jce - 0.0019142 * jce.powi(2) + jce.powi(3)/189474.0;

        // Mean anomaly of the sun (degrees)
        let x1 = 357.52772 + 35999.050340 * jce - 0.0001603 * jce.powi(2) - jce.powi(3)/300000.0;

        // Mean anomaly of the moon (degrees)
        let x2 = 134.96298 + 477198.867398 * jce + 0.0086972 * jce.powi(2) + jce.powi(3)/56250.0;

        // Moon's argument of latitude (degrees)
        let x3 = 93.27191 + 483202.017538 * jce - 0.0036825 * jce.powi(2) + jce.powi(3)/327270.0;

        // Moon's Longitude of ascending node (degrees)
        let x4 = 125.04452 - 1934.136261 * jce + 0.0020708 * jce.powi(2) + jce.powi(3)/450000.0;

        // Nutation in longitude and obliquity (degrees)
        let (delta_psi, delta_epsilon) = delta_psi_epsilon(x0, x1, x2, x3, x4, jce);
        
        // True obliquity of the ecliptic (degrees)
        let epsilon_deg = epsilon(jme, delta_epsilon);

        // Apparent sun longitude (degrees)
        let lambda_deg = longitude(theta_deg, delta_psi, r_au);

        // Apparent sidereal tim eat Greenwich (degrees)
        let nu_deg = nu(jd, jc, delta_psi, epsilon_deg);

        // Sun right ascention (degrees)
        let alpha_deg =  alpha(lambda_deg, epsilon_deg, beta_deg);

        // Geocentric sun declination (degrees)
        let delta_deg = delta(beta_deg, epsilon_deg, lambda_deg);

        // Observer local hour (degrees)
        let h_deg = h(nu_deg, obs_lon_deg, alpha_deg);

        // Equitorial parallax of the sun (degrees)
        let xi_deg = 8.794 / (3600.0 * r_au);

        // u-term (radians)
        let u_rad = (0.99664719 * obs_lat_deg.to_radians().tan()).atan();

        // Parallax in the sun right ascention (degrees)
        let delta_alpha_deg = delta_alpha(u_rad, obs_elev_m, obs_lat_deg, xi_deg, h_deg, delta_deg);

        // Topocentric sun right ascention (degrees)
        let alpha_prime_deg = alpha_deg + delta_alpha_deg;

        // Topocentric sun declination (degrees)
        let delta_prime_deg = delta_prime(u_rad, obs_elev_m, obs_lat_deg, delta_deg, xi_deg, delta_alpha_deg, h_deg);

        // Topocentric local hour (degrees)
        let h_prime_deg = h_deg - delta_alpha_deg;

        // Topocentric elevation angle (degrees)
        let e_deg = e(obs_lat_deg, delta_prime_deg, h_prime_deg, p_mbar, t_deg_c);

        // Topocentric zenith angle (degrees)
        let theta_small_deg = 90.0 - e_deg;

        // Topocentric astronomer's azimuth angle (degrees)
        let gamma_deg = gamma(h_prime_deg, obs_lat_deg, delta_prime_deg);

        // Topocentric azimuth angle (degrees)
        let phi_deg = (gamma_deg + 180.0) % 360.0;

        if DEBUG {
            println!("jd: {}", jd);
            println!("l: {}", l_deg);
            println!("b: {}", b_deg);
            println!("r: {}", r_au);
            println!("theta: {}", theta_deg);
            println!("beta: {}", beta_deg);
            println!("x0: {}", x0);
            println!("x1: {}", x1);
            println!("x2: {}", x2);
            println!("x3: {}", x3);
            println!("delta_phi: {}", delta_psi);
            println!("delta_epsilon: {}", delta_epsilon);
            println!("epsilon: {}", epsilon_deg);
            println!("lambda: {}", lambda_deg);
            println!("nu: {}", nu_deg);
            println!("alpha: {}", alpha_deg);
            println!("delta: {}", delta_deg);
            println!("H: {}", h_deg);
            println!("alpha_prime: {}", alpha_prime_deg);
            println!("H': {}", h_prime_deg);
            println!("delta_prime: {}", delta_prime_deg);
            println!("theta_small: {}", theta_small_deg);
            println!("phi: {}", phi_deg);
        }
    }
}

fn eq9(values: Vec<(f64, f64, f64)>, jme: f64) -> f64 {
    let mut s = 0.0;
    for v in values.iter() {
        s += v.0 * (v.1 + v.2 * jme).cos();
    }
    s
}

fn eq10(l0: f64, l1: f64, l2: f64, l3: f64, l4: f64, l5: f64, jme: f64) -> f64 {
    (l0 + l1*jme + l2*jme.powi(2) + l3*jme.powi(3) + l4*jme.powi(4) + l5*jme.powi(5)) / 10.0_f64.powi(8)
}

fn l_deg(jme: f64) -> f64 {
    let l0 = eq9(vec![
        (175347046.0, 0.0, 0.0),
        (3341656.0, 4.6692568, 6283.07585),
        (34894.0, 4.6261, 12566.1517),
        (3497.0, 2.7441, 5753.3849),
        (3418.0, 2.8289, 3.5231),
        (3136.0, 3.6277, 77713.7715),
        (2676.0, 4.4181, 7860.4194),
        (2343.0, 6.1352, 3930.2097),
        (1324.0, 0.7425, 11506.7698),
        (1273.0, 2.0371, 529.691),
        (1199.0, 1.1096, 1577.3435),
        (990.0, 5.233, 5884.927),
        (902.0, 2.045, 26.298),
        (857.0, 3.508, 398.149),
        (780.0, 1.179, 5223.694),
        (753.0, 2.533, 5507.553),
        (505.0, 4.583, 18849.228),
        (492.0, 4.205, 775.523),
        (357.0, 2.92, 0.067),
        (317.0, 5.849, 11790.629),
        (284.0, 1.899, 796.298),
        (271.0, 0.315, 10977.079),
        (243.0, 0.345, 5486.778),
        (206.0, 4.806, 2544.314),
        (205.0, 1.869, 5573.143),
        (202.0, 2.458, 6069.777),
        (156.0, 0.833, 213.299),
        (132.0, 3.411, 2942.463),
        (126.0, 1.083, 20.775),
        (115.0, 0.645, 0.98),
        (103.0, 0.636, 4694.003),
        (102.0, 0.976, 15720.839),
        (102.0, 4.267, 7.114),
        (99.0, 6.21, 2146.17),
        (98.0, 0.68, 155.42),
        (86.0, 5.98, 161000.69),
        (85.0, 1.3, 6275.96),
        (85.0, 3.67, 71430.7),
        (80.0, 1.81, 17260.15),
        (79.0, 3.04, 12036.46),
        (75.0, 1.76, 5088.63),
        (74.0, 3.5, 3154.69),
        (74.0, 4.68, 801.82),
        (70.0, 0.83, 9437.76),
        (62.0, 3.98, 8827.39),
        (61.0, 1.82, 7084.9),
        (57.0, 2.78, 6286.6),
        (56.0, 4.39, 14143.5),
        (56.0, 3.47, 6279.55),
        (52.0, 0.19, 12139.55),
        (52.0, 1.33, 1748.02),
        (51.0, 0.28, 5856.48),
        (49.0, 0.49, 1194.45),
        (41.0, 5.37, 8429.24),
        (41.0, 2.4, 19651.05),
        (39.0, 6.17, 10447.39),
        (37.0, 6.04, 10213.29),
        (37.0, 2.57, 1059.38),
        (36.0, 1.71, 2352.87),
        (36.0, 1.78, 6812.77),
        (33.0, 0.59, 17789.85),
        (30.0, 0.44, 83996.85),
        (30.0, 2.74, 1349.87),
        (25.0, 3.16, 4690.48)], jme);
    
    let l1 = eq9(vec![
        (628331966747.0, 0.0, 0.0),
        (206059.0, 2.678235, 6283.07585),
        (4303.0, 2.6351, 12566.1517),
        (425.0, 1.59, 3.523),
        (119.0, 5.796, 26.298),
        (109.0, 2.966, 1577.344),
        (93.0, 2.59, 18849.23),
        (72.0, 1.14, 529.69),
        (68.0, 1.87, 398.15),
        (67.0, 4.41, 5507.55),
        (59.0, 2.89, 5223.69),
        (56.0, 2.17, 155.42),
        (45.0, 0.4, 796.3),
        (36.0, 0.47, 775.52),
        (29.0, 2.65, 7.11),
        (21.0, 5.34, 0.98),
        (19.0, 1.85, 5486.78),
        (19.0, 4.97, 213.3),
        (17.0, 2.99, 6275.96),
        (16.0, 0.03, 2544.31),
        (16.0, 1.43, 2146.17),
        (15.0, 1.21, 10977.08),
        (12.0, 2.83, 1748.02),
        (12.0, 3.26, 5088.63),
        (12.0, 5.27, 1194.45),
        (12.0, 2.08, 4694.0),
        (11.0, 0.77, 553.57),
        (10.0, 1.3, 6286.6),
        (10.0, 4.24, 1349.87),
        (9.0, 2.7, 242.73),
        (9.0, 5.64, 951.72),
        (8.0, 5.3, 2352.87),
        (6.0, 2.65, 9437.76),
        (6.0, 4.67, 4690.48)], jme);
    
    let l2 = eq9(vec![
        (52919.0, 0.0, 0.0),
        (8720.0, 1.0721, 6283.0758),
        (309.0, 0.867, 12566.152),
        (27.0, 0.05, 3.52),
        (16.0, 5.19, 26.3),
        (16.0, 3.68, 155.42),
        (10.0, 0.76, 18849.23),
        (9.0, 2.06, 77713.77),
        (7.0, 0.83, 775.52),
        (5.0, 4.66, 1577.34),
        (4.0, 1.03, 7.11),
        (4.0, 3.44, 5573.14),
        (3.0, 5.14, 796.3),
        (3.0, 6.05, 5507.55),
        (3.0, 1.19, 242.73),
        (3.0, 6.12, 529.69),
        (3.0, 0.31, 398.15),
        (3.0, 2.28, 553.57),
        (2.0, 4.38, 5223.69),
        (2.0, 3.75, 0.98)], jme);

    let l3 = eq9(vec![
        (289.0, 5.844, 6283.076),
        (35.0, 0.0, 0.0),
        (17.0, 5.49, 12566.15),
        (3.0, 5.2, 155.42),
        (1.0, 4.72, 3.52),
        (1.0, 5.3, 18849.23),
        (1.0, 5.97, 242.73)], jme);

    let l4 = eq9(vec![
        (114.0, 3.142, 0.0),
        (8.0, 4.13, 6283.08),
        (1.0, 3.84, 12566.15)], jme);

    let l5 = 1.0 * (3.14 + 0.0 * jme).cos();

    let mut l_deg = eq10(l0, l1, l2, l3, l4, l5, jme).to_degrees() % 360.0;
    if l_deg < 0.0 {
        l_deg += 360.0;
    }
    l_deg
}

fn b_deg(jme: f64) -> f64 {
    let b0 = eq9(vec![
        (280.0, 3.199, 84334.662),
        (102.0, 5.422, 5507.553),
        (80.0, 3.88, 5223.69),
        (44.0, 3.7, 2352.87),
        (32.0, 4.0, 1577.34)], jme);
    let b1 = eq9(vec![(9.0, 3.9, 5507.55), (6.0, 1.73, 5223.69)], jme);
    eq10(b0, b1, 0.0, 0.0, 0.0, 0.0, jme).to_degrees()
}

fn r_au(jme: f64) -> f64 {
    let r0 = eq9(vec![
        (100013989.0, 0.0, 0.0),
        (1670700.0, 3.0984635, 6283.07585),
        (13956.0, 3.05525, 12566.1517),
        (3084.0, 5.1985, 77713.7715),
        (1628.0, 1.1739, 5753.3849),
        (1576.0, 2.8469, 7860.4194),
        (925.0, 5.453, 11506.77),
        (542.0, 4.564, 3930.21),
        (472.0, 3.661, 5884.927),
        (346.0, 0.964, 5507.553),
        (329.0, 5.9, 5223.694),
        (307.0, 0.299, 5573.143),
        (243.0, 4.273, 11790.629),
        (212.0, 5.847, 1577.344),
        (186.0, 5.022, 10977.079),
        (175.0, 3.012, 18849.228),
        (110.0, 5.055, 5486.778),
        (98.0, 0.89, 6069.78),
        (86.0, 5.69, 15720.84),
        (86.0, 1.27, 161000.69),
        (65.0, 0.27, 17260.15),
        (63.0, 0.92, 529.69),
        (57.0, 2.01, 83996.85),
        (56.0, 5.24, 71430.7),
        (49.0, 3.25, 2544.31),
        (47.0, 2.58, 775.52),
        (45.0, 5.54, 9437.76),
        (43.0, 6.01, 6275.96),
        (39.0, 5.36, 4694.0),
        (38.0, 2.39, 8827.39),
        (37.0, 0.83, 19651.05),
        (37.0, 4.9, 12139.55),
        (36.0, 1.67, 12036.46),
        (35.0, 1.84, 2942.46),
        (33.0, 0.24, 7084.9),
        (32.0, 0.18, 5088.63),
        (32.0, 1.78, 398.15),
        (28.0, 1.21, 6286.6),
        (28.0, 1.9, 6279.55),
        (26.0, 4.59, 10447.39)], jme);
    
    let r1 = eq9(vec![
        (103019.0, 1.10749, 6283.07585),
        (1721.0, 1.0644, 12566.1517),
        (702.0, 3.142, 0.0),
        (32.0, 1.02, 18849.23),
        (31.0, 2.84, 5507.55),
        (25.0, 1.32, 5223.69),
        (18.0, 1.42, 1577.34),
        (10.0, 5.91, 10977.08),
        (9.0, 1.42, 6275.96),
        (9.0, 0.27, 5486.78)], jme);
    
    let r2 = eq9(vec![
        (4359.0, 5.7846, 6283.0758),
        (124.0, 5.579, 12566.152),
        (12.0, 3.14, 0.0),
        (9.0, 3.63, 77713.77),
        (6.0, 1.87, 5573.14),
        (3.0, 5.47, 18849.23)], jme);
    let r3 = eq9(vec![(145.0, 4.273, 6283.076), (7.0, 3.92, 12566.15)], jme);
    let r4 = eq9(vec![(4.0, 2.56, 6283.08)], jme);
    eq10(r0, r1, r2, r3, r4, 0.0, jme)
}

fn delta_psi_epsilon(x0: f64, x1: f64, x2: f64, x3: f64, x4: f64, jce: f64) -> (f64, f64) {
    let values = vec![
        (0.0, 0.0, 0.0, 0.0, 1.0, -171996.0, -174.2, 92025.0, 8.9),
        (-2.0, 0.0, 0.0, 2.0, 2.0, -13187.0, -1.6, 5736.0, -3.1),
        (0.0, 0.0, 0.0, 2.0, 2.0, -2274.0, -0.2, 977.0, -0.5),
        (0.0, 0.0, 0.0, 0.0, 2.0, 2062.0, 0.2, -895.0, 0.5),
        (0.0, 1.0, 0.0, 0.0, 0.0, 1426.0, -3.4, 54.0, -0.1),
        (0.0, 0.0, 1.0, 0.0, 0.0, 712.0, 0.1, -7.0, 0.0),
        (-2.0, 1.0, 0.0, 2.0, 2.0, -517.0, 1.2, 224.0, -0.6),
        (0.0, 0.0, 0.0, 2.0, 1.0, -386.0, -0.4, 200.0, 0.0),
        (0.0, 0.0, 1.0, 2.0, 2.0, -301.0, 0.0, 129.0, -0.1),
        (-2.0, -1.0, 0.0, 2.0, 2.0, 217.0, -0.5, -95.0, 0.3),
        (-2.0, 0.0, 1.0, 0.0, 0.0, -158.0, 0.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0, 2.0, 1.0, 129.0, 0.1, -70.0, 0.0),
        (0.0, 0.0, -1.0, 2.0, 2.0, 123.0, 0.0, -53.0, 0.0),
        (2.0, 0.0, 0.0, 0.0, 0.0, 63.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0, 1.0, 63.0, 0.1, -33.0, 0.0),
        (2.0, 0.0, -1.0, 2.0, 2.0, -59.0, 0.0, 26.0, 0.0),
        (0.0, 0.0, -1.0, 0.0, 1.0, -58.0, -0.1, 32.0, 0.0),
        (0.0, 0.0, 1.0, 2.0, 1.0, -51.0, 0.0, 27.0, 0.0),
        (-2.0, 0.0, 2.0, 0.0, 0.0, 48.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, -2.0, 2.0, 1.0, 46.0, 0.0, -24.0, 0.0),
        (2.0, 0.0, 0.0, 2.0, 2.0, -38.0, 0.0, 16.0, 0.0),
        (0.0, 0.0, 2.0, 2.0, 2.0, -31.0, 0.0, 13.0, 0.0),
        (0.0, 0.0, 2.0, 0.0, 0.0, 29.0, 0.0, 0.0, 0.0),
        (-2.0, 0.0, 1.0, 2.0, 2.0, 29.0, 0.0, -12.0, 0.0),
        (0.0, 0.0, 0.0, 2.0, 0.0, 26.0, 0.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0, 2.0, 0.0, -22.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, -1.0, 2.0, 1.0, 21.0, 0.0, -10.0, 0.0),
        (0.0, 2.0, 0.0, 0.0, 0.0, 17.0, -0.1, 0.0, 0.0),
        (2.0, 0.0, -1.0, 0.0, 1.0, 16.0, 0.0, -8.0, 0.0),
        (-2.0, 2.0, 0.0, 2.0, 2.0, -16.0, 0.1, 7.0, 0.0),
        (0.0, 1.0, 0.0, 0.0, 1.0, -15.0, 0.0, 9.0, 0.0),
        (-2.0, 0.0, 1.0, 0.0, 1.0, -13.0, 0.0, 7.0, 0.0),
        (0.0, -1.0, 0.0, 0.0, 1.0, -12.0, 0.0, 6.0, 0.0),
        (0.0, 0.0, 2.0, -2.0, 0.0, 11.0, 0.0, 0.0, 0.0),
        (2.0, 0.0, -1.0, 2.0, 1.0, -10.0, 0.0, 5.0, 0.0),
        (2.0, 0.0, 1.0, 2.0, 2.0, -8.0, 0.0, 3.0, 0.0),
        (0.0, 1.0, 0.0, 2.0, 2.0, 7.0, 0.0, -3.0, 0.0),
        (-2.0, 1.0, 1.0, 0.0, 0.0, -7.0, 0.0, 0.0, 0.0),
        (0.0, -1.0, 0.0, 2.0, 2.0, -7.0, 0.0, 3.0, 0.0),
        (2.0, 0.0, 0.0, 2.0, 1.0, -7.0, 0.0, 3.0, 0.0),
        (2.0, 0.0, 1.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0),
        (-2.0, 0.0, 2.0, 2.0, 2.0, 6.0, 0.0, -3.0, 0.0),
        (-2.0, 0.0, 1.0, 2.0, 1.0, 6.0, 0.0, -3.0, 0.0),
        (2.0, 0.0, -2.0, 0.0, 1.0, -6.0, 0.0, 3.0, 0.0),
        (2.0, 0.0, 0.0, 0.0, 1.0, -6.0, 0.0, 3.0, 0.0),
        (0.0, -1.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0),
        (-2.0, -1.0, 0.0, 2.0, 1.0, -5.0, 0.0, 3.0, 0.0),
        (-2.0, 0.0, 0.0, 0.0, 1.0, -5.0, 0.0, 3.0, 0.0),
        (0.0, 0.0, 2.0, 2.0, 1.0, -5.0, 0.0, 3.0, 0.0),
        (-2.0, 0.0, 2.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0),
        (-2.0, 1.0, 0.0, 2.0, 1.0, 4.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, -2.0, 0.0, 4.0, 0.0, 0.0, 0.0),
        (-1.0, 0.0, 1.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0),
        (-2.0, 1.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, -2.0, 2.0, 2.0, -3.0, 0.0, 0.0, 0.0),
        (-1.0, -1.0, 1.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 1.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0),
        (0.0, -1.0, 1.0, 2.0, 2.0, -3.0, 0.0, 0.0, 0.0),
        (2.0, -1.0, -1.0, 2.0, 2.0, -3.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 3.0, 2.0, 2.0, -3.0, 0.0, 0.0, 0.0),
        (2.0, -1.0, 0.0, 2.0, 2.0, -3.0, 0.0, 0.0, 0.0)];
    let mut delta_psi = 0.0;
    let mut delta_epsilon = 0.0;
    for v in values.iter() {
        let s = x0 * v.0 + x1 * v.1 + x2 * v.2 + x3 * v.3 + x4 * v.4;
        delta_psi += (v.5 + v.6 * jce) * s.to_radians().sin();
        delta_epsilon += (v.7 + v.8 * jce) * s.to_radians().cos();
    }
    return (delta_psi / 36000000.0, delta_epsilon / 36000000.0);
}

// Eq. 24 & 25 (degrees)
fn epsilon(jme: f64, delta_epsilon: f64) -> f64 {
    let u = jme / 10.0;
    let epsilon0 = 84381.448 - 4680.93*u - 1.55*u.powi(2) + 1999.25*u.powi(3) -
    51.38*u.powi(4) - 249.67*u.powi(5) - 39.05*u.powi(6) + 7.12*u.powi(7) +
    27.87*u.powi(8) + 5.79*u.powi(9) + 2.45*u.powi(10);
    epsilon0/3600.0 + delta_epsilon
}

// Eq. 26 & 27 (degrees)
fn longitude(theta_deg: f64, delta_psi: f64, r: f64) -> f64{
    let delta_tau = -20.4898 / (3600.0 * r);
    theta_deg + delta_psi + delta_tau
}

// Eq. 28 & 29 (degrees)
fn nu(jd: f64, jc: f64, delta_psi: f64, epsilon_deg: f64) -> f64 {
    let mut nu0 = (280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*jc.powi(2) - jc.powi(3)/38710000.0) % 360.0;
    if nu0 < 0.0 {
        nu0 += 360.0;
    }
    nu0 + delta_psi*epsilon_deg.to_radians().cos()
}

// Eq. 30 (degrees)
fn alpha(lambda_deg: f64, epsilon_deg: f64, beta_deg: f64) -> f64 {
    let mut alpha = (lambda_deg.to_radians().sin()*epsilon_deg.to_radians().cos() - 
    beta_deg.to_radians().tan()*epsilon_deg.to_radians().sin()).atan2(lambda_deg.to_radians().cos());
    alpha = alpha.to_degrees() % 360.0;
    if alpha < 0.0 {
        alpha += 360.0;
    }
    alpha
}

// Eq. 31 (degrees)
fn delta(beta_deg: f64, epsilon_deg: f64, lambda_deg: f64) -> f64 {
    (beta_deg.to_radians().sin()*epsilon_deg.to_radians().cos() + 
        beta_deg.to_radians().cos()*epsilon_deg.to_radians().sin()*lambda_deg.to_radians().sin()).asin().to_degrees()
}

// Eq. 32 (degrees)
fn h(nu_deg: f64, obs_lon_deg: f64, alpha_deg: f64) -> f64 {
    let mut h = (nu_deg + obs_lon_deg - alpha_deg) % 360.0;
    if h < 0.0 {
        h += 360.0
    }
    h
}

// Eq. 37 (degrees)
fn delta_alpha(u_rad: f64, obs_elev_m: f64, obs_lat_deg: f64, xi_deg: f64, h_deg: f64, delta_deg: f64) -> f64 {
    let x = u_rad.cos() + (obs_elev_m / 6378140.0) * obs_lat_deg.to_radians().cos();
    (-x * xi_deg.to_radians().sin() * h_deg.to_radians().sin()).atan2(
        delta_deg.to_radians().cos() - x * xi_deg.to_radians().sin() * h_deg.to_radians().cos()).to_degrees()
}

// Eq. 39 (degrees)
fn delta_prime(u_rad: f64, obs_elev_m: f64, obs_lat_deg: f64, delta_deg: f64, xi_deg: f64, delta_alpha_deg: f64, h_deg: f64) -> f64 {
    let x = u_rad.cos() + (obs_elev_m / 6378140.0) * obs_lat_deg.to_radians().cos();
    let y = 0.99664719 * u_rad.sin() + (obs_elev_m / 6378140.0) * obs_lat_deg.to_radians().sin();
    ((delta_deg.to_radians().sin() - y * xi_deg.to_radians().sin()) * delta_alpha_deg.to_radians().cos()).atan2(
        delta_deg.to_radians().cos() - x * xi_deg.to_radians().sin() * h_deg.to_radians().cos()).to_degrees()
}

fn e(obs_lat_deg: f64, delta_prime_deg: f64, h_prime_deg: f64, p_mbar: f64, t_deg_c: f64) -> f64 {
    let e0 = (obs_lat_deg.to_radians().sin() * delta_prime_deg.to_radians().sin() + 
        obs_lat_deg.to_radians().cos() * delta_prime_deg.to_radians().cos() * h_prime_deg.to_radians().cos()).asin().to_degrees();
    let delta_e = (p_mbar / 1010.0) * (283.0 / (273.0 + t_deg_c)) * (1.02 / (60.0 * (e0 + 10.3 / (e0 + 5.11)).to_radians().tan()));
    e0 + delta_e
}

fn gamma(h_prime_deg: f64, obs_lat_deg: f64, delta_prime_deg: f64) -> f64 {
   let mut lambda = h_prime_deg.to_radians().sin().atan2(
        h_prime_deg.to_radians().cos() * obs_lat_deg.to_radians().sin() - delta_prime_deg.to_radians().tan() * obs_lat_deg.to_radians().cos()).to_degrees();
    if lambda < 0.0 {
        lambda += 360.0;
    }
    lambda
}