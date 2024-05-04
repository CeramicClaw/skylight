use crate::*;

pub const DEBUG: bool = false;

pub struct Moon;

impl Moon {
    // Get a characteristic Sun at a given DateTime, observer lat/lon, elevation, pressure, and temperature
    pub fn new_moon(placetime: &PlaceTime) -> CelestialObject {
        // Julian dates
        let jd = julian_day(&placetime.date);
        let jde = julian_ephemeris_day(jd, placetime.date.delta_t);
        let jc = julian_century(jd);
        let jce = julian_ephemeris_century(jde);
        let jme = jce / 10.0;

        // l, r, and b terms
        let d = moon_mean_elongation(jce);
        let m = sun_mean_anomaly(jce);
        let m_prime = moon_mean_anomaly(jce);
        let f = moon_arg_of_latitude(jce);
        let l_prime = moon_mean_longitude(jce);
        let l_term = l_term(d, m, m_prime, f, jce);
        let r_term = r_term(d, m, m_prime, f, jce);
        let b_term = b_term(d, m, m_prime, f, jce);

        // Delta l and Delta b
        let a1 = a1(jce); // Eq. 18 (degrees)
        let a2 = a2(jce); // Eq. 19 (degrees)
        let a3 = a3(jce); // Eq. 20 (degrees)
        let delta_l = delta_l(a1, a2, l_prime, f);
        let delta_b = delta_b(a1, a3, l_prime, m_prime, f);

        // Moon longitude and latitude (degrees)
        let moon_lon = moon_longitude(l_prime, l_term, delta_l);
        let beta = moon_latitutde(b_term, delta_b);

        // Moon distance from center of Earth
        let delta = moon_distance(r_term);

        // Horizontal parallax
        let pi = moon_horizontal_parallax(delta); // Radians

        // Nutation in Longitude and Obliquity
        let x0 = x0(jce);
        let x1 = x1(jce);
        let x2 = x2(jce);
        let x3 = x3(jce);
        let x4 = x4(jce);
        let delta_psi = delta_psi(jce, x0, x1, x2, x3, x4);
        let delta_epsilon = delta_epsilon(jce, x0, x1, x2, x3, x4);

        // True Obliquity of the Ecliptic
        let epsilon = epsilon(jme, delta_epsilon);
        
        // Apparent Moon Longitude
        let lambda = moon_lon + delta_psi; // Eq. 38

        // Apparent Sidereal Time at Greenwich
        let eta0 = eta0(jd, jc);
        let eta = eta(eta0, delta_psi, epsilon);

        // Moon's Right Ascension
        let alpha = alpha(lambda, epsilon, beta);

        // Moon's Geocentric Declination
        let delta_small = delta_small(beta, epsilon, lambda);

        // Observer Local Hour Angle
        let h = obs_local_angle(eta, placetime.obs_lon_deg, alpha);

        // Moon's Topocentric Right Ascention
        let u = u(placetime.obs_lat_deg); // Radians
        let x = x(u, placetime.obs_lat_deg, placetime.obs_elev_m);
        let y = y(u, placetime.obs_lat_deg, placetime.obs_elev_m);
        let delta_alpha = delta_alpha(x, pi, h, delta_small);
        let alpha_prime = alpha + delta_alpha; // Eq. 48

        // Moon's Topocentric Declination
        let delta_small_prime = delta_small_prime(delta_small, x, y, pi, delta_alpha, h);
        
        // Topocentric Local Hour (degrees)
        let h_prime = h - delta_alpha; // Eq. 50

        // Topocentric Elevation Angle Without Correction
        let e0 = e0(placetime.obs_lat_deg, delta_small_prime, h_prime);

        // Atmospheric Refraction correction
        let delta_e = delta_e(placetime.p_mbar, placetime.t_deg_c, e0);
        let e = e0 + delta_e;
        let theta_m = 90.0 - e; // Eq. 54
        let gamma = gamma(h_prime, placetime.obs_lat_deg, delta_small_prime);
        
        // Topocentric Azimuth
        let phi_m = (gamma + 180.0) % 360.0;

        // Radius of Moon's Disk
        let r_m = r_m(e, pi, delta);

        if DEBUG {
            println!("JD: {}\nJCE: {}\nD: {}\nM: {}\nM': {}\nF: {}\nL': {}", jd, jce, d, m, m_prime, f, l_prime);
            println!("l_term: {}", l_term);
            println!("r_term: {}", r_term);
            println!("b_term: {}", b_term);
            println!("a1: {}\na2: {}\na3: {}", a1, a2, a3);
            println!("l + delta_l: {}", l_term + delta_l);
            println!("b + delta_b: {}", b_term + delta_b);
            println!("lambda prime: {}", moon_lon);
            println!("Beta: {}", beta);
            println!("Delta: {}", delta);
            println!("Pi: {}", pi);
            println!("delta_psi: {}", delta_psi);
            println!("delta_epsilon: {}", delta_epsilon);
            println!("epsilon: {}", epsilon);
            println!("lambda: {}", lambda);
            println!("eta: {}", eta);
            println!("H: {}", h);
            println!("delta_alpha: {}", delta_alpha);
            println!("h_prime {}", h_prime);
            println!("delta_e {}", delta_e);
            println!("e {}", e);
            println!("theta_m {}", theta_m);
            println!("gamma {}", gamma);
            println!("phi_m {}", phi_m);
            println!("===========");
            println!("Horizontal parallax: {}", pi.to_degrees());
            println!("Geocentric Right Ascention: {}", alpha);
            println!("Geocentric Declination: {}", delta_small);
            println!("Topcentric Right Ascention: {}", alpha_prime);
            println!("Topcentric Declination: {}", delta_small_prime);
            println!("Moon Apparent Radius: {}", r_m);
        }        

        return CelestialObject {
            date: placetime.date,
            horiz_parallax_deg: pi.to_degrees(),
            geo_r_asc_deg: alpha,
            geo_dec_deg: delta_small,
            topo_r_asc_deg: alpha_prime,
            topo_dec_deg: delta_small_prime,
            semidiameter_deg: r_m,
            zenith_deg: theta_m,
            azimuth_deg: phi_m,
            celestial: Celestial::MOON,
        };
    }
}

/// Eq. 9 (degrees)
fn moon_mean_longitude(t: f64) -> f64 {
    let mut l_prime = 218.3164477 + 481267.88123421*t - 0.0015786*t.powi(2) + t.powi(3)/538841.0 
        - t.powi(4)/65194000.0;
    l_prime = l_prime % 360.0;
    if l_prime < 0.0 {
        l_prime += 360.0;
    }
    l_prime
}

/// Eq. 10 (degrees)
fn moon_mean_elongation(t: f64) -> f64 {
    let mut d = 297.8501921 + 445267.1114034*t - 0.0018819*t.powi(2) + t.powi(3)/545868.0 
    - t.powi(4)/113065000.0;
    d = d % 360.0;
    if d < 0.0 {
        d += 360.0;
    }
    d
}

/// Eq. 11 (degrees)
fn sun_mean_anomaly(t: f64) -> f64 {
    357.5291092 + 35999.0502909*t - 0.0001536*t.powi(2) + t.powi(3)/24490000.0
}

/// Eq. 12 (degrees)
fn moon_mean_anomaly(t: f64) -> f64 {
    134.9633964 + 477198.8675055*t + 0.0087414*t.powi(2) + t.powi(3)/69699.0 - t.powi(4)/14712000.0
}

// Eq. 13 (degrees)
fn moon_arg_of_latitude(t: f64) -> f64 {
    93.2720950 + 483202.0175233*t - 0.0036539*t.powi(2) - t.powi(3)/3526000.0 + t.powi(4)/863310000.0
}

// Eq. 14 (0.000001 degrees)
fn l_term(d: f64, m: f64, m_prime: f64, f: f64, t: f64) -> f64 {
    let values: Vec<(i32, i32, i32, i32, i32)> = vec![
        (0, 0, 1, 0, 6288774),
        (2, 0, -1, 0, 1274027),
        (2, 0, 0, 0, 658314),
        (0, 0, 2, 0, 213618),
        (0, 1, 0, 0, -185116),
        (0, 0, 0, 2, -114332),
        (2, 0, -2, 0, 58793),
        (2, -1, -1, 0, 57066),
        (2, 0, 1, 0, 53322),
        (2, -1, 0, 0, 45758),
        (0, 1, -1, 0, -40923),
        (1, 0, 0, 0, -34720),
        (0, 1, 1, 0, -30383),
        (2, 0, 0, -2, 15327),
        (0, 0, 1, 2, -12528),
        (0, 0, 1, -2, 10980),
        (4, 0, -1, 0, 10675),
        (0, 0, 3, 0, 10034),
        (4, 0, -2, 0, 8548),
        (2, 1, -1, 0, -7888),
        (2, 1, 0, 0, -6766),
        (1, 0, -1, 0, -5163),
        (1, 1, 0, 0, 4987),
        (2, -1, 1, 0, 4036),
        (2, 0, 2, 0, 3994),
        (4, 0, 0, 0, 3861),
        (2, 0, -3, 0, 3665),
        (0, 1, -2, 0, -2689),
        (2, 0, -1, 2, -2602),
        (2, -1, -2, 0, 2390),
        (1, 0, 1, 0, -2348),
        (2, -2, 0, 0, 2236),
        (0, 1, 2, 0, -2120),
        (0, 2, 0, 0, -2069),
        (2, -2, -1, 0, 2048),
        (2, 0, 1, -2, -1773),
        (2, 0, 0, 2, -1595),
        (4, -1, -1, 0, 1215),
        (0, 0, 2, 2, -1110),
        (3, 0, -1, 0, -892),
        (2, 1, 1, 0, -810),
        (4, -1, -2, 0, 759),
        (0, 2, -1, 0, -713),
        (2, 2, -1, 0, -700),
        (2, 1, -2, 0, 691),
        (2, -1, 0, -2, 596),
        (4, 0, 1, 0, 549),
        (0, 0, 4, 0, 537),
        (4, -1, 0, 0, 520),
        (1, 0, -2, 0, -487),
        (2, 1, 0, -2, -399),
        (0, 0, 2, -2, -381),
        (1, 1, 1, 0, 351),
        (3, 0, -2, 0, -340),
        (4, 0, -3, 0, 330),
        (2, -1, 2, 0, 327),
        (0, 2, 1, 0, -323),
        (1, 1, -1, 0, 299),
        (2, 0, 3, 0, 294)];

    let e = 1.0 - 0.002516 * t - 0.0000074 * t * t;
    let mut l = 0.0;
    let mut l_i;
    for v in values.iter() {
        l_i = v.4 as f64;
        if v.1.abs() == 1 {
            l_i *= e;
        }
        if v.1.abs() == 2 {
            l_i *= e * e;
        }
        l += l_i * (v.0 as f64 * d + v.1 as f64 * m + v.2 as f64 * m_prime + v.3 as f64 * f).to_radians().sin();
    }
    l
}

/// Eq. 16 (0.001 km)
fn r_term(d: f64, m: f64, m_prime: f64, f: f64, t: f64) -> f64 {
    let values: Vec<(i32, i32, i32, i32, i32)> = vec![
        (0, 0, 1, 0, -20905355),
        (2, 0, -1, 0, -3699111),
        (2, 0, 0, 0, -2955968),
        (0, 0, 2, 0, -569925),
        (0, 1, 0, 0, 48888),
        (0, 0, 0, 2, -3149),
        (2, 0, -2, 0, 246158),
        (2, -1, -1, 0, -152138),
        (2, 0, 1, 0, -170733),
        (2, -1, 0, 0, -204586),
        (0, 1, -1, 0, -129620),
        (1, 0, 0, 0, 108743),
        (0, 1, 1, 0, 104755),
        (2, 0, 0, -2, 10321),
        (0, 0, 1, -2, 79661),
        (4, 0, -1, 0, -34782),
        (0, 0, 3, 0, -23210),
        (4, 0, -2, 0, -21636),
        (2, 1, -1, 0, 24208),
        (2, 1, 0, 0, 30824),
        (1, 0, -1, 0, -8379),
        (1, 1, 0, 0, -16675),
        (2, -1, 1, 0, -12831),
        (2, 0, 2, 0, -10445),
        (4, 0, 0, 0, -11650),
        (2, 0, -3, 0, 14403),
        (0, 1, -2, 0, -7003),
        (2, -1, -2, 0, 10056),
        (1, 0, 1, 0, 6322),
        (2, -2, 0, 0, -9884),
        (0, 1, 2, 0, 5751),
        (2, -2, -1, 0, -4950),
        (2, 0, 1, -2, 4130),
        (4, -1, -1, 0, -3958),
        (3, 0, -1, 0, 3258),
        (2, 1, 1, 0, 2616),
        (4, -1, -2, 0, -1897),
        (0, 2, -1, 0, -2117),
        (2, 2, -1, 0, 2354),
        (4, 0, 1, 0, -1423),
        (0, 0, 4, 0, -1117),
        (4, -1, 0, 0, -1571),
        (1, 0, -2, 0, -1739),
        (0, 0, 2, -2, -4421),
        (0, 2, 1, 0, 1165),
        (2, 0, -1, -2, 8752)];

    let e = 1.0 - 0.002516 * t - 0.0000074 * t * t;
    let mut r = 0.0;
    let mut r_i;
    for v in values.iter() {
        r_i = v.4 as f64;
        if v.1.abs() == 1 {
            r_i *= e;
        }
        if v.1.abs() == 2 {
            r_i *= e * e;
        }
        r += r_i * (v.0 as f64 * d + v.1 as f64 * m + v.2 as f64 * m_prime + v.3 as f64 * f).to_radians().cos();
    }
    r
}

fn b_term(d: f64, m: f64, m_prime: f64, f: f64, t: f64) -> f64 {
    let values: Vec<(i32, i32, i32, i32, i32)> = vec![
        (0, 0, 0, 1, 5128122),
        (0, 0, 1, 1, 280602),
        (0, 0, 1, -1, 277693),
        (2, 0, 0, -1, 173237),
        (2, 0, -1, 1, 55413),
        (2, 0, -1, -1, 46271),
        (2, 0, 0, 1, 32573),
        (0, 0, 2, 1, 17198),
        (2, 0, 1, -1, 9266),
        (0, 0, 2, -1, 8822),
        (2, -1, 0, -1, 8216),
        (2, 0, -2, -1, 4324),
        (2, 0, 1, 1, 4200),
        (2, 1, 0, -1, -3359),
        (2, -1, -1, 1, 2463),
        (2, -1, 0, 1, 2211),
        (2, -1, -1, -1, 2065),
        (0, 1, -1, -1, -1870),
        (4, 0, -1, -1, 1828),
        (0, 1, 0, 1, -1794),
        (0, 0, 0, 3, -1749),
        (0, 1, -1, 1, -1565),
        (1, 0, 0, 1, -1491),
        (0, 1, 1, 1, -1475),
        (0, 1, 1, -1, -1410),
        (0, 1, 0, -1, -1344),
        (1, 0, 0, -1, -1335),
        (0, 0, 3, 1, 1107),
        (4, 0, 0, -1, 1021),
        (4, 0, -1, 1, 833),
        (0, 0, 1, -3, 777),
        (4, 0, -2, 1, 671),
        (2, 0, 0, -3, 607),
        (2, 0, 2, -1, 596),
        (2, -1, 1, -1, 491),
        (2, 0, -2, 1, -451),
        (0, 0, 3, -1, 439),
        (2, 0, 2, 1, 422),
        (2, 0, -3, -1, 421),
        (2, 1, -1, 1, -366),
        (2, 1, 0, 1, -351),
        (4, 0, 0, 1, 331),
        (2, -1, 1, 1, 315),
        (2, -2, 0, -1, 302),
        (0, 0, 1, 3, -283),
        (2, 1, 1, -1, -229),
        (1, 1, 0, -1, 223),
        (1, 1, 0, 1, 223),
        (0, 1, -2, -1, -220),
        (2, 1, -1, -1, -220),
        (1, 0, 1, 1, -185),
        (2, -1, -2, -1, 181),
        (0, 1, 2, 1, -177),
        (4, 0, -2, -1, 176),
        (4, -1, -1, -1, 166),
        (1, 0, 1, -1, -164),
        (4, 0, 1, -1, 132),
        (1, 0, -1, -1, -119),
        (4, -1, 0, -1, 115),
        (2, -2, 0, 1, 107)];

    let e = 1.0 - 0.002516 * t - 0.0000074 * t * t;
    let mut b = 0.0;
    let mut b_i;
    for v in values.iter() {
        b_i = v.4 as f64;
        if v.1.abs() == 1 {
            b_i *= e;
        }
        if v.1.abs() == 2 {
            b_i *= e * e;
        }
        b += b_i * (v.0 as f64 * d + v.1 as f64 * m + v.2 as f64 * m_prime + v.3 as f64 * f).to_radians().sin();
    }
    b
}

// Eq. 18 (degrees)
fn a1(t: f64) -> f64 {
    119.75 + 131.849 * t
}

// Eq. 19 (degrees)
fn a2(t: f64) -> f64 {
    53.09 + 479264.29 * t
}

// Eq. 20 (degrees)
fn a3(t: f64) -> f64 {
    313.45 + 481266.484 * t
}

// Eq. 21 (degrees)
fn delta_l(a1: f64, a2: f64, l_prime: f64, f: f64) -> f64{
    3958.0 * a1.to_radians().sin() + 1962.0 * (l_prime - f).to_radians().sin() + 
    318.0 * a2.to_radians().sin()
}

// Eq. 22 (degrees)
fn delta_b(a1: f64, a3: f64, l_prime: f64, m_prime: f64, f: f64) -> f64 {
    -2235.0 * l_prime.to_radians().sin() + 382.0 * a3.to_radians().sin() + 
    175.0 * (a1 - f).to_radians().sin() + 175.0 * (a1 + f).to_radians().sin() + 
    127.0 * (l_prime - m_prime).to_radians().sin() - 115.0 * (l_prime + m_prime).to_radians().sin()
}

// Eq. 23 (degrees)
fn moon_longitude(l_prime: f64, l_term: f64, delta_l: f64) -> f64 {
    let mut moon_lon = (l_prime + (l_term + delta_l) / 1000000.0) % 360.0;
    moon_lon = moon_lon % 360.0;
    if moon_lon < 0.0 {
        moon_lon += 360.0;
    }
    moon_lon
}

// Eq. 24 (degrees)
fn moon_latitutde(b_term: f64, delta_b: f64) -> f64 {
    let mut moon_lat = (b_term + delta_b) / 1000000.0;
    moon_lat = moon_lat % 360.0;
    if moon_lat < 0.0 {
        moon_lat += 360.0;
    }
    moon_lat
}

// Eq. 25 (km)
fn moon_distance(r: f64) -> f64 {
    385000.56 + r / 1000.0
}

// Eq. 26 (radians)
fn moon_horizontal_parallax(delta: f64) -> f64 {
    (6378.14 / delta).asin()
}

// Eq. 27 (degrees)
fn x0(jce: f64) -> f64 {
    297.85036 + 445267.111480 * jce - 0.0019142 * jce.powi(2) + jce.powi(3)/189474.0
}

// Eq. 28 (degrees)
fn x1(jce: f64) -> f64 {
    357.52772 + 35999.050340 * jce - 0.0001603 * jce.powi(2) + jce.powi(3)/300000.0
}

// Eq. 29 (degrees)
fn x2(jce: f64) -> f64 {
    134.96298 + 477198.867398 * jce + 0.0086972 * jce.powi(2) + jce.powi(3)/56250.0
}

// Eq. 30 (degrees)
fn x3(jce: f64) -> f64 {
    93.27191 + 483202.017538 * jce - 0.0036825 * jce.powi(2) + jce.powi(3)/327270.0
}

// Eq. 31 (degrees)
fn x4(jce: f64) -> f64 {
    125.04452 - 1934.136261 * jce + 0.0020708 * jce.powi(2) + jce.powi(3)/450000.0
}

// Eq. 32 & 34 (degrees)
fn delta_psi(jce: f64, x0: f64, x1: f64, x2: f64, x3: f64, x4: f64) -> f64 {
    // Y0, Y1, Y2, Y3, Y4, a, B
    let values: Vec<(i32, i32, i32, i32, i32, i32, f64)> = vec![
        (0, 0, 0, 0, 1, -171996, -174.2),
        (-2, 0, 0, 2, 2, -13187, -1.6),
        (0, 0, 0, 2, 2, -2274, -0.2),
        (0, 0, 0, 0, 2, 2062, 0.2),
        (0, 1, 0, 0, 0, 1426, -3.4),
        (0, 0, 1, 0, 0, 712, 0.1),
        (-2, 1, 0, 2, 2, -517, 1.2),
        (0, 0, 0, 2, 1, -386, -0.4),
        (0, 0, 1, 2, 2, -301, 0.0),
        (-2, -1, 0, 2, 2, 217, -0.5),
        (-2, 0, 1, 0, 0, -158, 0.0),
        (-2, 0, 0, 2, 1, 129, 0.1),
        (0, 0, -1, 2, 2, 123, 0.0),
        (2, 0, 0, 0, 0, 63, 0.0),
        (0, 0, 1, 0, 1, 63, 0.1),
        (2, 0, -1, 2, 2, -59, 0.0),
        (0, 0, -1, 0, 1, -58, -0.1),
        (0, 0, 1, 2, 1, -51, 0.0),
        (-2, 0, 2, 0, 0, 48, 0.0),
        (0, 0, -2, 2, 1, 46, 0.0),
        (2, 0, 0, 2, 2, -38, 0.0),
        (0, 0, 2, 2, 2, -31, 0.0),
        (0, 0, 2, 0, 0, 29, 0.0),
        (-2, 0, 1, 2, 2, 29, 0.0),
        (0, 0, 0, 2, 0, 26, 0.0),
        (-2, 0, 0, 2, 0, -22, 0.0),
        (0, 0, -1, 2, 1, 21, 0.0),
        (0, 2, 0, 0, 0, 17, -0.1),
        (2, 0, -1, 0, 1, 16, 0.0),
        (-2, 2, 0, 2, 2, -16, 0.1),
        (0, 1, 0, 0, 1, -15, 0.0),
        (-2, 0, 1, 0, 1, -13, 0.0),
        (0, -1, 0, 0, 1, -12, 0.0),
        (0, 0, 2, -2, 0, 11, 0.0),
        (2, 0, -1, 2, 1, -10, 0.0),
        (2, 0, 1, 2, 2, -8, 0.0),
        (0, 1, 0, 2, 2, 7, 0.0),
        (-2, 1, 1, 0, 0, -7, 0.0),
        (0, -1, 0, 2, 2, -7, 0.0),
        (2, 0, 0, 2, 1, -7, 0.0),
        (2, 0, 1, 0, 0, 6, 0.0),
        (-2, 0, 2, 2, 2, 6, 0.0),
        (-2, 0, 1, 2, 1, 6, 0.0),
        (2, 0, -2, 0, 1, -6, 0.0),
        (2, 0, 0, 0, 1, -6, 0.0),
        (0, -1, 1, 0, 0, 5, 0.0),
        (-2, -1, 0, 2, 1, -5, 0.0),
        (-2, 0, 0, 0, 1, -5, 0.0),
        (0, 0, 2, 2, 1, -5, 0.0),
        (-2, 0, 2, 0, 1, 4, 0.0),
        (-2, 1, 0, 2, 1, 4, 0.0),
        (0, 0, 1, -2, 0, 4, 0.0),
        (-1, 0, 1, 0, 0, -4, 0.0),
        (-2, 1, 0, 0, 0, -4, 0.0),
        (1, 0, 0, 0, 0, -4, 0.0),
        (0, 0, 1, 2, 0, 3, 0.0),
        (0, 0, -2, 2, 2, -3, 0.0),
        (-1, -1, 1, 0, 0, -3, 0.0),
        (0, 1, 1, 0, 0, -3, 0.0),
        (0, -1, 1, 2, 2, -3, 0.0),
        (2, -1, -1, 2, 2, -3, 0.0),
        (0, 0, 3, 2, 2, -3, 0.0),
        (2, -1, 0, 2, 2, -3, 0.0)];
    
    let mut delta_psi = 0.0;
    for v in values.iter() {
        delta_psi += (v.5 as f64 + v.6 * jce) * 
            (x0 * v.0 as f64 + x1 * v.1 as f64 + x2 * v.2 as f64 + x3 * v.3 as f64 + x4 * v.4 as f64).to_radians().sin();
    }
    delta_psi / 36000000.0
}

// Eq. 33 & 35 (degrees)
fn delta_epsilon(jce: f64, x0: f64, x1: f64, x2: f64, x3: f64, x4: f64) -> f64 {
    // Y0, Y1, Y2, Y3, Y4, c, d
    let values: Vec<(i32, i32, i32, i32, i32, i32, f64)> = vec![
        (0, 0, 0, 0, 1, 92025, 8.9),
        (-2, 0, 0, 2, 2, 5736, -3.1),
        (0, 0, 0, 2, 2, 977, -0.5),
        (0, 0, 0, 0, 2, -895, 0.5),
        (0, 1, 0, 0, 0, 54, -0.1),
        (0, 0, 1, 0, 0, -7, 0.0),
        (-2, 1, 0, 2, 2, 224, -0.6),
        (0, 0, 0, 2, 1, 200, 0.0),
        (0, 0, 1, 2, 2, 129, -0.1),
        (-2, -1, 0, 2, 2, -95, 0.3),
        (-2, 0, 0, 2, 1, -70, 0.0),
        (0, 0, -1, 2, 2, -53, 0.0),
        (0, 0, 1, 0, 1, -33, 0.0),
        (2, 0, -1, 2, 2, 26, 0.0),
        (0, 0, -1, 0, 1, 32, 0.0),
        (0, 0, 1, 2, 1, 27, 0.0),
        (0, 0, -2, 2, 1, -24, 0.0),
        (2, 0, 0, 2, 2, 16, 0.0),
        (0, 0, 2, 2, 2, 13, 0.0),
        (-2, 0, 1, 2, 2, -12, 0.0),
        (0, 0, -1, 2, 1, -10, 0.0),
        (2, 0, -1, 0, 1, -8, 0.0),
        (-2, 2, 0, 2, 2, 7, 0.0),
        (0, 1, 0, 0, 1, 9, 0.0),
        (-2, 0, 1, 0, 1, 7, 0.0),
        (0, -1, 0, 0, 1, 6, 0.0),
        (2, 0, -1, 2, 1, 5, 0.0),
        (2, 0, 1, 2, 2, 3, 0.0),
        (0, 1, 0, 2, 2, -3, 0.0),
        (0, -1, 0, 2, 2, 3, 0.0),
        (2, 0, 0, 2, 1, 3, 0.0),
        (-2, 0, 2, 2, 2, -3, 0.0),
        (-2, 0, 1, 2, 1, -3, 0.0),
        (2, 0, -2, 0, 1, 3, 0.0),
        (2, 0, 0, 0, 1, 3, 0.0),
        (-2, -1, 0, 2, 1, 3, 0.0),
        (-2, 0, 0, 0, 1, 3, 0.0),
        (0, 0, 2, 2, 1, 3, 0.0)];
    
    let mut delta_epsilon = 0.0;
    for v in values.iter() {
        delta_epsilon += (v.5 as f64 + v.6 * jce) * 
            (x0 * v.0 as f64 + x1 * v.1 as f64 + x2 * v.2 as f64 + x3 * v.3 as f64 + x4 * v.4 as f64).to_radians().cos();
    }
    delta_epsilon / 36000000.0
}

// Eq. 36 & 37 (degrees)
fn epsilon(jme: f64, delta_epsilon: f64) -> f64 {
    let u = jme / 10.0;
    let epsilon0 = 84381.448 - 4680.93*u - 1.55*u.powi(2) + 1999.25*u.powi(3) -
    51.38*u.powi(4) - 249.67*u.powi(5) - 39.05*u.powi(6) + 7.12*u.powi(7) +
    27.87*u.powi(8) + 5.79*u.powi(9) + 2.45*u.powi(10);
    epsilon0/3600.0 + delta_epsilon
}

// Eq. 39
fn eta0(jd: f64, jc: f64) -> f64 {
    280.46061837 + 360.98564736629*(jd - 2451545.0) + 0.000387933*jc.powi(2) - jc.powi(3)/38710000.0
}

// Eq. 40
fn eta(eta0: f64, delta_psi: f64, epsilon: f64) -> f64 {
    let mut eta = eta0 + delta_psi*epsilon.to_radians().cos();
    eta = eta % 360.0;
    if eta < 0.0 {
        eta += 360.0;
    }
    eta
}

// Eq 41 (degrees)
fn alpha(lambda: f64, epsilon: f64, beta: f64) -> f64 {
    let mut alpha = (lambda.to_radians().sin()*epsilon.to_radians().cos() - 
        beta.to_radians().tan()*epsilon.to_radians().sin()).atan2(lambda.to_radians().cos());
    alpha = alpha.to_degrees() % 360.0;
    if alpha < 0.0 {
        alpha += 360.0;
    }
    alpha
}

// Eq. 42 (degrees)
fn delta_small(beta: f64, epsilon: f64, lambda: f64) -> f64 {
    (beta.to_radians().sin()*epsilon.to_radians().cos() + 
        beta.to_radians().cos()*epsilon.to_radians().sin()*lambda.to_radians().sin()).asin().to_degrees()
}

// Eq. 43
fn obs_local_angle(eta: f64, obs_lon: f64, alpha: f64) -> f64 {
    let mut h = (eta + obs_lon - alpha) % 360.0;
    if h < 0.0 {
        h += 360.0;
    }
    h
}

// Eq. 44 (Radians)
fn u(obs_lat: f64) -> f64 {
    (0.99664719 * obs_lat.to_radians().tan()).atan()
}

// Eq. 45
fn x(u: f64, obs_lat: f64, obs_elev: f64) -> f64 {
    u.cos() + (obs_elev * obs_lat.to_radians().cos())/6378140.0
}

// Eq. 46
fn y(u: f64, obs_lat: f64, obs_elev: f64) -> f64 {
    0.99664719*u.sin() + (obs_elev * obs_lat.to_radians().sin())/6378140.0
}

// Eq. 47 (degrees)
fn delta_alpha(x: f64, pi: f64, h: f64, delta_small: f64) -> f64 {
    ((-x*pi.sin()*h.to_radians().sin()).atan2(delta_small.to_radians().cos() - 
        x*pi.sin()*h.to_radians().cos())).to_degrees()
}

// Eq. 49 (degrees)
fn delta_small_prime(delta_small: f64, x: f64, y: f64, pi: f64, delta_alpha: f64, h: f64) -> f64 {
    (((delta_small.to_radians().sin() - y*pi.sin())*delta_alpha.to_radians().cos()).atan2(
        delta_small.to_radians().cos() - x*pi.sin()*h.to_radians().cos())).to_degrees()
}

// Eq. 51 (degrees)
fn e0(obs_lat: f64, delta_small_prime: f64, h_prime: f64) -> f64 {
    (obs_lat.to_radians().sin() * delta_small_prime.to_radians().sin() +
        obs_lat.to_radians().cos() * delta_small_prime.to_radians().cos()*h_prime.to_radians().cos()).asin().to_degrees()
}

// Eq. 52 (degrees)
// p: pressure in millibars
// t: temperature in degrees C
fn delta_e(p: f64, t: f64, e0: f64) -> f64 {
    (p/1010.0) * (283.0/(273.0 + t)) * 1.02 / (60.0 * ((e0 + 10.3 / (e0 + 5.11))).to_radians().tan())
}

// Eq. 55 (degrees)
fn gamma(h_prime: f64, obs_lat: f64, delta_small_prime: f64) -> f64 {
    let mut lambda = h_prime.to_radians().sin().atan2(
        h_prime.to_radians().cos()*obs_lat.to_radians().sin() - 
        delta_small_prime.to_radians().tan()*obs_lat.to_radians().cos()).to_degrees();
    lambda = lambda % 360.0;
    if lambda < 0.0 {
        lambda += 360.0;
    }
    lambda
}

// Eq. 59 (degrees)
fn r_m(e: f64, pi: f64, delta: f64) -> f64 {
    358473400.0 * (1.0 + e.to_radians().sin() * pi.sin()) / (3600.0 * delta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Month;

    #[test]
    fn test() { // Values from 1998 Astronomical Almanac
        let mut pt = PlaceTime::new(&new_day_delta(1998, Month::JANUARY, 5, 0.0), 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut moon = Moon::new_moon(&pt);
        assert!((0.98754 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((0.6908 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::JANUARY, 10, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((0.962467 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((17.649167 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::FEBRUARY, 1, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((1.000856 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((-0.692 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::FEBRUARY, 5, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((0.96860278 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((14.98194 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::MAY, 17, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((0.95081389 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((-17.3805 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::AUGUST, 16, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((0.976783 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((16.50472 - moon.geo_dec_deg).abs() < 0.003);
        pt.change_time(&new_day_delta(1998, Month::AUGUST, 20, 0.0));
        moon = Moon::new_moon(&pt);
        assert!((0.944694 - moon.horiz_parallax_deg).abs() < 0.001);
        assert!((17.151 - moon.geo_dec_deg).abs() < 0.003);
    }
}