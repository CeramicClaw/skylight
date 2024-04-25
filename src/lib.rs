use core::fmt;

pub mod moon;
pub mod sun;

pub const DELTA_T: f64 = 69.0; // Default delta_t value as of January 2024

#[derive(PartialEq, Debug)]
pub enum Sign {
    POSITIVE,
    NEGATIVE
}

pub enum Month {
    JANUARY,
    FEBRUARY,
    MARCH,
    APRIL,
    MAY,
    JUNE,
    JULY,
    AUGUST,
    SEPTEMBER,
    OCTOBER,
    NOVEMBER,
    DECEMBER,
}

pub struct DateTime {
    pub year: i32,
    pub month: Month,
    pub day: u32,
    pub hour: u32,
    pub minute: u32,
    pub second: f64,
    pub delta_t: f64,
}

#[derive(PartialEq, Debug)]
pub struct DMS {
    pub degree: u32,
    pub minute: u32,
    pub second: f64,
    pub sign: Sign,
}

#[derive(PartialEq, Debug)]
pub struct HMS {
    pub hour: u32,
    pub minute: u32,
    pub second: f64,
}

impl HMS {
    pub fn decimal(&self) -> f64 {
        self.hour as f64 * (360.0 / 24.0) + 
        self.minute as f64 * (360.0 / (24.0 * 60.0)) +
        self.second * (360.0 / (24.0 * 60.0 * 60.0))
    }
}

impl DMS {
    pub fn decimal(&self) -> f64 {
        println!("{}", self);
        let d = self.degree as f64 + self.minute as f64 / 60.0 + self.second / 3600.0;
        if self.sign == Sign::POSITIVE {
            return d;
        } else {
            return d * -1.0;
        }
    }
}

impl Month {
    pub fn value(&self) -> u8 {
        match self {
            Month::JANUARY => 1,
            Month::FEBRUARY => 2,
            Month::MARCH => 3,
            Month::APRIL => 4,
            Month::MAY => 5,
            Month::JUNE => 6,
            Month::JULY => 7,
            Month::AUGUST => 8,
            Month::SEPTEMBER => 9,
            Month::OCTOBER => 10,
            Month::NOVEMBER => 11,
            Month::DECEMBER => 12,
        }
    }
    pub fn text(&self) -> &str {
        match self {
            Month::JANUARY => "January",
            Month::FEBRUARY => "February",
            Month::MARCH => "March",
            Month::APRIL => "April",
            Month::MAY => "May",
            Month::JUNE => "June",
            Month::JULY => "July",
            Month::AUGUST => "August",
            Month::SEPTEMBER => "September",
            Month::OCTOBER => "October",
            Month::NOVEMBER => "November",
            Month::DECEMBER => "December",
        }
    }
}

impl DateTime {
    /// Get the fractional decimal day value
    pub fn decimal_day(&self) -> f64 {
        self.day as f64 + 
        self.hour as f64 / 24.0 + // Hours in a day
        self.minute as f64 / 1440.0 + // Minutes in a day
        self.second / 86400.0 // Seconds in a day
    }
}

impl std::fmt::Display for DateTime {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}, {} ({}:{}:{})",
        self.month.text(), 
        self.day,
        self.year,
        format!("{:0>2}", self.hour),
        format!("{:0>2}", self.minute),
        format!("{:0>2}", self.second),
        )
    }
}

impl std::fmt::Display for DMS {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}\u{00B0} {}\u{2032} {:.2}\u{2033}",
        if self.sign == Sign::NEGATIVE {
            "-"
        } else {
            ""
        },
        self.degree,
        self.minute,
        self.second,
        )
    }
}

impl std::fmt::Display for HMS {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\u{02B0} {}\u{1D50} {:.2}\u{02E2}",
        self.hour,
        self.minute,
        self.second,
        )
    }
}

pub fn dms(d: i32, m: i32, s: f64) -> DMS {
    let mut sign = Sign::POSITIVE;
    if (d as f64) < 0.0 || (m as f64) < 0.0 || s < 0.0 {
        sign = Sign::NEGATIVE;
    }
    DMS{degree: d.abs() as u32, minute: m.abs() as u32, second: s.abs(), sign}
}

pub fn hms(hour: u32, minute: u32, second: f64) -> HMS {
    HMS{hour, minute, second}
}

pub fn deg2dms(d: f64) -> DMS {
    let mut sign = Sign::POSITIVE;
    if d < 0.0 {
        sign  = Sign::NEGATIVE;
    }
    let degree = d.abs().trunc() as u32;
    let minute = ((d.abs() - degree as f64) * 60.0).trunc() as u32;
    let second = (d.abs() - degree as f64 - minute as f64 / 60.0) * 3600.0;
    DMS{degree, minute, second, sign}
}

pub fn deg2hms(degree: f64) -> HMS {
    if degree < 0.0 {
        panic!("Unable to convert a negative value ({})to HMS format!", degree);
    }
    let hour = (degree * 24.0 / 360.0).trunc() as u32;
    let minute = ((degree - hour as f64 * 360.0 / 24.0) * (24.0 * 60.0) / 360.0).trunc() as u32;
    let second = (degree - hour as f64 * 360.0 / 24.0 - minute as f64 * 360.0 / (24.0 * 60.0)) * (24.0 * 60.0 * 60.0) / 360.0;
    HMS{hour, minute, second}
}

pub fn new_day(year: i32, month: Month, day: u32) -> DateTime {
    DateTime {
        year,
        month,
        day,
        hour: 0,
        minute: 0,
        second: 0.0,
        delta_t: DELTA_T,
    }
}

pub fn new_time(year: i32, month: Month, day: u32, hour: u32, minute: u32, second: f64) -> DateTime {
    DateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
        delta_t: DELTA_T,
    }
}

pub fn new_time_t(year: i32, month: Month, day: u32, hour: u32, minute: u32, second: f64, delta_t: f64) -> DateTime {
    DateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
        delta_t,
    }
}

/// Eq. 4
pub fn julian_day(date: &DateTime) -> f64 {
    let mut y = date.year as f64;
    let mut m = date.month.value() as f64;
    if date.month.value() == 1 || date.month.value() == 2 {
        y = y - 1.0;
        m = m + 12.0;
    }
    let d = date.decimal_day();
    let jd = (365.25 * (y + 4716.0)).trunc() + (30.6001 * (m + 1.0)).trunc() + d - 1524.5;
    let mut b = 0.0;
    if jd > 2299160.0 {
        let a = (y / 100.0).trunc();
        b = 2.0 - a + (a / 4.0).trunc();
    }
    return jd + b;
}

/// Eq. 5
pub fn julian_ephemeris_day(jd: f64, delta_t: f64) -> f64 {
    jd + delta_t / 86400.0
}

/// Eq. 6
pub fn julian_century(jd: f64) -> f64 {
    (jd - 2451545.0) / 36525.0
}

/// Eq. 7
pub fn julian_ephemeris_century(jde: f64) -> f64 {
    (jde - 2451545.0) / 36525.0
}

/// Eq. 8
pub fn julian_ephemeris_millenium(jce: f64) -> f64 {
    jce / 10.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Month;

    #[test]
    fn test() { // Table A4.1 test values
        assert_eq!(2451545.0, julian_day(&new_time(2000, Month::JANUARY, 1, 12, 0, 0.0)));
        assert_eq!(2451179.5, julian_day(&new_time(1999, Month::JANUARY, 1, 0, 0, 0.0)));
        assert_eq!(2446822.5, julian_day(&new_time(1987, Month::JANUARY, 27, 0, 0, 0.0)));
        assert_eq!(2446966.0, julian_day(&new_time(1987, Month::JUNE, 19, 12, 0, 0.0)));
        assert_eq!(2447187.5, julian_day(&new_time(1988, Month::JANUARY, 27, 0, 0, 0.0)));
        assert_eq!(2447332.0, julian_day(&new_time(1988, Month::JUNE, 19, 12, 0, 0.0)));
        assert_eq!(2415020.5, julian_day(&new_time(1900, Month::JANUARY, 1, 0, 0, 0.0)));
        assert_eq!(2305447.5, julian_day(&new_time(1600, Month::JANUARY, 1, 0, 0, 0.0)));
        assert_eq!(2305812.5, julian_day(&new_time(1600, Month::DECEMBER, 31, 0, 0, 0.0)));
        assert_eq!(2026871.8, julian_day(&new_time(837, Month::APRIL, 10, 7, 12, 0.0)));
        assert_eq!(1676496.5, julian_day(&new_time(-123, Month::DECEMBER, 31, 0, 0, 0.0)));
        assert_eq!(1676497.5, julian_day(&new_time(-122, Month::JANUARY, 1, 0, 0, 0.0)));
        assert_eq!(1356001.0, julian_day(&new_time(-1000, Month::JULY, 12, 12, 0, 0.0)));
        assert_eq!(1355866.5, julian_day(&new_time(-1000, Month::FEBRUARY, 29, 0, 0, 0.0)));
        assert_eq!(1355671.4, julian_day(&new_time(-1001, Month::AUGUST, 17, 21, 36, 0.0)));
        assert_eq!(0.0, julian_day(&new_time(-4712, Month::JANUARY, 1, 12, 0, 0.0)));
        // Sun Example
        assert_eq!(2452930.312847222, julian_day(&new_time_t(2003, Month::OCTOBER, 17, 19, 30, 30.0, 67.0)));
        // DMS
        assert_eq!(0.0, dms(0, 0, 0.0).decimal());
        assert_eq!(90.0, dms(90, 0, 0.0).decimal());
        assert_eq!(-170.0, dms(-170, 0, 0.0).decimal());
        assert!(0.0000001 > 179.99972222 - dms(179, 59, 59.0).decimal());
        assert_eq!(-30.255555555555556, dms(-30, 15, 20.0).decimal());
        assert!(0.00000001 > (dms(-30, 15, 20.0).decimal() - deg2dms(dms(-30, 15, 20.0).decimal()).decimal()).abs());
        // HMS
        assert_eq!(0.0, hms(0, 0, 0.0).decimal());
        assert_eq!(45.0, hms(3, 0, 0.0).decimal());
        assert_eq!(5.1875, hms(0, 20, 45.0).decimal());
        assert!(0.00000001 > (hms(10, 15, 20.0).decimal() - deg2hms(hms(10, 15, 20.0).decimal()).decimal()).abs());
    }
}