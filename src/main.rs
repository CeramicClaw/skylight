use skylight::*;
use skylight::moon::*;
use skylight::sun::*;

fn main() {
    let t = new_time(2024, Month::APRIL, 8, 18, 31, 4.0);
    let pt = PlaceTime::new(&t, 46.062, -68.74, 0.0, 820.0, 11.0);
    let moon = Moon::new_moon(&pt);
    let sun = Sun::new_day(&pt);
    println!("{}", sun);
    println!("{}", moon);
    println!("Is eclipse? {}", is_solar_eclipse(sun, moon));
}
