use core::fmt;

pub mod moon;

pub const DELTA_T: f64 = 69.0; // Default delta_t value as of January 2024

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

pub struct DateTime {
    pub year: i32,
    pub month: Month,
    pub day: u32,
    pub hour: u32,
    pub minute: u32,
    pub second: u32,
    pub delta_t: f64,
}

impl DateTime {
    /// Get the fractional decimal day value
    pub fn decimal_day(&self) -> f64 {
        self.day as f64 + 
        self.hour as f64 / 24.0 + // Hours in a day
        self.minute as f64 / 1440.0 + // Minutes in a day
        self.second as f64 / 3600.0 // Seconds in a day
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

pub fn new_day(year: i32, month: Month, day: u32) -> DateTime {
    DateTime {
        year,
        month,
        day,
        hour: 0,
        minute: 0,
        second: 0,
        delta_t: DELTA_T,
    }
}

pub fn new_time(year: i32, month: Month, day: u32, hour: u32, minute: u32, second: u32) -> DateTime {
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

pub fn new_time_t(year: i32, month: Month, day: u32, hour: u32, minute: u32, second: u32, delta_t: f64) -> DateTime {
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
        assert_eq!(2451545.0, julian_day(&new_time(2000, Month::JANUARY, 1, 12, 0, 0)));
        assert_eq!(2451179.5, julian_day(&new_time(1999, Month::JANUARY, 1, 0, 0, 0)));
        assert_eq!(2446822.5, julian_day(&new_time(1987, Month::JANUARY, 27, 0, 0, 0)));
        assert_eq!(2446966.0, julian_day(&new_time(1987, Month::JUNE, 19, 12, 0, 0)));
        assert_eq!(2447187.5, julian_day(&new_time(1988, Month::JANUARY, 27, 0, 0, 0)));
        assert_eq!(2447332.0, julian_day(&new_time(1988, Month::JUNE, 19, 12, 0, 0)));
        assert_eq!(2415020.5, julian_day(&new_time(1900, Month::JANUARY, 1, 0, 0, 0)));
        assert_eq!(2305447.5, julian_day(&new_time(1600, Month::JANUARY, 1, 0, 0, 0)));
        assert_eq!(2305812.5, julian_day(&new_time(1600, Month::DECEMBER, 31, 0, 0, 0)));
        assert_eq!(2026871.8, julian_day(&new_time(837, Month::APRIL, 10, 7, 12, 0)));
        assert_eq!(1676496.5, julian_day(&new_time(-123, Month::DECEMBER, 31, 0, 0, 0)));
        assert_eq!(1676497.5, julian_day(&new_time(-122, Month::JANUARY, 1, 0, 0, 0)));
        assert_eq!(1356001.0, julian_day(&new_time(-1000, Month::JULY, 12, 12, 0, 0)));
        assert_eq!(1355866.5, julian_day(&new_time(-1000, Month::FEBRUARY, 29, 0, 0, 0)));
        assert_eq!(1355671.4, julian_day(&new_time(-1001, Month::AUGUST, 17, 21, 36, 0)));
        assert_eq!(0.0, julian_day(&new_time(-4712, Month::JANUARY, 1, 12, 0, 0)));
    }
}