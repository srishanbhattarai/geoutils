//! Geoutils is a evolving crate to provide several geological computations and utilities.
//! Most computations are based off methods on the [Location](struct.Location.html) struct. A simple example to
//! get distance between two points:
//! ```
//! extern crate geoutils;
//!
//! use geoutils::Location;
//!
//! let berlin = Location::new(52.518611, 13.408056);
//! let moscow = Location::new(55.751667, 37.617778);
//! let distance = berlin.distance_to(&moscow).unwrap();
//!
//! println!("Distance = {}", distance);
//! ```

#![deny(missing_docs)]
mod formula;

/// Location defines a point using it's latitude and longitude.
pub struct Location(f64, f64);

impl Location {
    /// Create a new Location with it's degree values of latitude and longitude.
    pub fn new<T: Into<f64>>(lat: T, lon: T) -> Self {
        Location(lat.into(), lon.into())
    }

    /// Get the latitude.
    pub fn latitude(&self) -> f64 {
        self.0
    }

    /// Get the longitude.
    pub fn longitude(&self) -> f64 {
        self.1
    }

    /// Find the distance from itself to another point. Internally uses Vincenty's inverse formula.
    /// For better performance and lesser accuracy, consider [haversine_distance_to](struct.Location.html#method.haversine_distance_to).
    pub fn distance_to(&self, to: &Location) -> Result<f64, String> {
        match formula::vincenty_inverse(self, to, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance),
            Err(e) => Err(e),
        }
    }

    /// Find the distance from itself to another point using Haversine formula.
    /// This is usually computationally less intensive than [distance_to](struct.Location.html#method.distance_to) but
    /// is generally not as accurate.
    pub fn haversine_distance_to(&self, to: &Location) -> f64 {
        formula::haversine_distance_to(self, to)
    }

    /// Check if the point is within a fixed radius of another point.
    pub fn is_in_circle<T: Into<f64>>(&self, center: &Location, radius: T) -> Result<bool, String> {
        match formula::vincenty_inverse(self, center, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance < radius.into()),
            Err(e) => Err(e),
        }
    }

    /// Find the center of given locations.
    pub fn center(coords: Vec<&Location>) -> Location {
        formula::center_of_coords(coords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_get_distance() {
        let l1 = Location::new(27.740068, 85.337576);
        let l2 = Location::new(27.740286, 85.337059);

        match l1.distance_to(&l2) {
            Ok(distance) => {
                assert_eq!(distance, 56.409);
            }
            Err(e) => panic!("Failed: {:?}", e),
        }
    }

    #[test]
    fn test_get_distance_haversine() {
        let l1 = Location::new(27.740068, 85.337576);
        let l2 = Location::new(27.740286, 85.337059);

        let distance = l1.haversine_distance_to(&l2);
        assert_eq!(distance, 56.36);
    }

    #[test]
    fn test_get_center() {
        let l1 = Location::new(52.518611, 13.408056);
        let l2 = Location::new(55.751667, 37.617778);

        let l = Location::center(vec![&l1, &l2]);

        assert_eq!(l.0, 54.743683);
        assert_eq!(l.1, 25.033239);
    }
}
