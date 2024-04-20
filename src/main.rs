use skylight::*;
//use skylight::moon::*;
use skylight::sun::*;

fn main() {
    let t = new_time_t(2003, Month::OCTOBER, 17, 19, 30, 30, 67.0);
    //let mut moon = Moon::new_moon();
    //moon.set(t, 24.61167, 143.36167);
    let mut sun = Sun::new_day();
    sun.set(t, 39.742476, -105.1786, 1830.14, 820.0, 11.0);
    println!("{}", sun);
    // Page ~200 has the moon positions
}
