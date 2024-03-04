use skylight::*;
use skylight::moon::*;

fn main() {
    let mut t = new_time(2009, Month::JULY, 22, 1, 33, 0);
    //let mut t = new_time(1992, Month::APRIL, 12, 0, 0, 0);
    t.delta_t = 66.4;
    println!("{}", Moon::new_moon(t, 24.61167, 143.36167));
}
