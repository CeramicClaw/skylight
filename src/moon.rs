use crate::*;

pub struct Moon {
    date: DateTime,
}

impl std::fmt::Display for Moon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}",
        self.date,
        )
    }
}

impl Moon {
    /// Get a characteristic moon at a given DateTime
    pub fn new_moon(date: DateTime) -> Moon {
        let mut moon = Moon::default();
        moon.date = date;
        // Calculate Julian dates
        let jd = julian_day(&moon.date);
        let jde = julian_ephemeris_day(jd, moon.date.delta_t);
        //let jc = julian_century(jd);
        let jce = julian_ephemeris_century(jde);
        //let jce = -0.077221081451;

        // Calculate l, r, and b terms
        let d = moon_mean_elongation(jce);
        let m = sun_mean_anomaly(jce);
        let m_prime = moon_mean_anomaly(jce);
        let f = moon_arg_of_latitude(jce);
        let l_prime = moon_mean_longitude(jce);
        println!("JD: {}\nJCE: {}\nD: {}\nM: {}\nM': {}\nF: {}\nL': {}", jd, jce, d, m, m_prime, f, l_prime);
        let l_term = l_term(d, m, m_prime, f, jce);
        println!("l_term: {}", l_term);
        let r_term = r_term(d, m, m_prime, f, jce);
        println!("r_term: {}", r_term);
        let b_term = b_term(d, m, m_prime, f, jce);
        println!("b_term: {}", b_term);

        // Calculate delta l and delta b
        let a1 = a1(jce); // Eq. 18 (degrees)
        let a2 = a2(jce); // Eq. 19 (degrees)
        let a3 = a3(jce); // Eq. 20 (degrees)
        println!("a1: {}\na2: {}\na3: {}", a1, a2, a3);
        let delta_l = delta_l(a1, a2, l_prime, f);
        println!("l + delta_l: {}", l_term + delta_l);
        let delta_b = delta_b(a1, a3, l_prime, m_prime, f);
        println!("b + delta_b: {}", b_term + delta_b);

        // Calculate moon longitude and latitude (degrees)
        let moon_lon = moon_longitude(l_prime, l_term, delta_l);
        println!("lambda prime: {}", moon_lon);
        let moon_lat = moon_latitutde(b_term, delta_b);
        println!("Beta: {}", moon_lat);

        let delta = moon_distance(r_term);
        println!("Delta: {}", delta);

        let pi = moon_horizontal_parallax(delta);
        println!("Pi: {}", pi);

        moon
    }

    fn default() -> Moon {
        Moon {
            date: new_day(2000, Month::JANUARY, 1),
        }
    }
}

/// Eq. 9 (degrees)
fn moon_mean_longitude(t: f64) -> f64 {
    let mut l_prime = 218.3164477 + 481267.88123421 * t - 0.0015786 * t * t + (t * t * t) / 538841.0 - (t * t * t * t) / 65194000.0;
    l_prime = l_prime % 360.0;
    if l_prime < 0.0 {
        l_prime += 360.0;
    }
    l_prime
}

/// Eq. 10 (degrees)
fn moon_mean_elongation(t: f64) -> f64 {
    let mut d = 297.8501921 + 445267.1114034 * t - 0.0018819 * t * t + (t * t * t) / 545868.0 - (t * t * t * t) /113065000.0;
    d = d % 360.0;
    if d < 0.0 {
        d += 360.0;
    }
    d
}

/// Eq. 11 (degrees)
fn sun_mean_anomaly(t: f64) -> f64 {
    let m = 357.5291092 + 35999.0502909 * t - 0.0001536 * t * t + (t *t * t) / 24490000.0;
    m
}

/// Eq. 12 (degrees)
fn moon_mean_anomaly(t: f64) -> f64 {
    let m_prime = 134.9633964 + 477198.8675055 * t + 0.0087414 * t * t + (t * t * t) / 69699.0 - (t * t * t * t) /14712000.0;
    m_prime
}

// Eq. 13 (degrees)
fn moon_arg_of_latitude(t: f64) -> f64 {
    let f = 93.2720950 + 483202.0175233 * t - 0.0036539 * t * t - (t * t * t) / 3526000.0 + (t * t * t * t) / 863310000.0;
    f
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
    let a1 = 119.75 + 131.849 * t;
    a1
}

// Eq. 19 (degrees)
fn a2(t: f64) -> f64 {
    let a2 = 53.09 + 479264.29 * t;
    a2
}

// Eq. 20 (degrees)
fn a3(t: f64) -> f64 {
    let a3 = 313.45 + 481266.484 * t;
    a3
}

// Eq. 21 (degrees)
fn delta_l(a1: f64, a2: f64, l_prime: f64, f: f64) -> f64{
    let delta_l = 3958.0 * a1.to_radians().sin() + 1962.0 * (l_prime - f).to_radians().sin() + 
    318.0 * a2.to_radians().sin();
    delta_l
}

// Eq. 22 (degrees)
fn delta_b(a1: f64, a3: f64, l_prime: f64, m_prime: f64, f: f64) -> f64 {
    let delta_b = -2235.0 * l_prime.to_radians().sin() + 382.0 * a3.to_radians().sin() + 
    175.0 * (a1 - f).to_radians().sin() + 175.0 * (a1 + f).to_radians().sin() + 
    127.0 * (l_prime - m_prime).to_radians().sin() - 115.0 * (l_prime + m_prime).to_radians().sin();
    delta_b
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
    let dist = 385000.56 + r / 1000.0;
    dist
}

// Eq. 26 (radians)
fn moon_horizontal_parallax(delta: f64) -> f64 {
    let pi = (6378.14 / delta).asin();
    pi
}