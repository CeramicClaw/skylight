use skylight::*;
use skylight::moon::*;

fn main() {
    //let mut t = new_time(2009, Month::JULY, 22, 1, 33, 0);
    let mut t = new_time(1998, Month::FEBRUARY, 6, 0, 0, 0);
    t.delta_t = 66.4;
    //t.delta_t = 0.0;
    let mut moon = Moon::new_moon();
    moon.set(t, 24.61167, 143.36167);
    println!("{}", moon);
    // Page ~200 has the moon positions
    //println!("{}", Moon::new_moon(t, 0.0, 0.0));
}
