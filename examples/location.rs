extern crate geoutils;

use geoutils::Location;

fn main() {
    // Create two locations
    let berlin = Location::new(52.518611, 13.408056);
    let moscow = Location::new(55.751667, 37.617778);

    // Get distance between two locations
    let distance = berlin.distance_to(&moscow).unwrap();
    println!("Distance between berlin and moscow: {} meters", distance);

    // Get the center of a number of locations
    let center = Location::center(vec![&berlin, &moscow]);
    println!(
        "Center of berlin and moscow: Lat={}, Lon={}",
        center.0, center.1
    );
}
