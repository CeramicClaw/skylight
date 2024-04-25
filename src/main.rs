use skylight::*;
use skylight::moon::*;
use skylight::sun::*;

fn main() {
    //let t = new_time_t(2003, Month::OCTOBER, 17, 19, 30, 30.0, 67.0);
    let t1 = new_time(2024, Month::APRIL, 8, 18, 31, 4.0);
    let t2 = new_time(2024, Month::APRIL, 8, 18, 31, 4.0);
    //let d = hms(10, 59, 2.1);
    //println!("d: {}", deg2hms(d.decimal()));
    let mut moon = Moon::new_moon();
    moon.set(t1, 46.062, -68.74);
    let mut sun = Sun::new_day();
    sun.set(t2, 46.062, -68.74, 0.0, 820.0, 11.0);
    println!("{}", sun);
    println!("{}", moon);
    // Page ~200 has the moon positions
}
