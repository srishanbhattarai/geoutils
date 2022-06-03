//! Geoutils is a evolving crate to provide several geological computations and utilities.
//! Most computations are based off methods on the [Location](struct.Location.html) struct.
//!
//! Find the full API reference at [docs.rs](http://docs.rs/geoutils).
//!
//! # Examples
//!
//! * Get the distance between two points using [Vincenty's Inverse Formula](https://en.wikipedia.org/wiki/Vincenty%27s_formulae).
//! ```rust
//! use geoutils::Location;
//!
//! let berlin = Location::new(52.518611, 13.408056);
//! let moscow = Location::new(55.751667, 37.617778);
//! let distance = berlin.distance_to(&moscow).unwrap();
//!
//! println!("Distance = {}", distance.meters());
//! ```
//!
//! * Get the distance between two points using the [Haversine Formula](https://en.wikipedia.org/wiki/Haversine_formula).
//! ```rust
//! use geoutils::Location;
//!
//! let berlin = Location::new(52.518611, 13.408056);
//! let moscow = Location::new(55.751667, 37.617778);
//! let distance = berlin.haversine_distance_to(&moscow);
//!
//! println!("Distance = {}", distance.meters());
//! ```
//!
//! * Get the center of a list of coordinates.
//! ```rust
//! use geoutils::Location;
//!
//! let berlin = Location::new(52.518611, 13.408056);
//! let moscow = Location::new(55.751667, 37.617778);
//! let center = Location::center(&vec![&berlin, &moscow]);
//!
//! println!("Center {}, {}", center.latitude(), center.longitude());
//! ```
//!
//! * Check if a point falls in a certain radius of another point.
//! ```rust
//! use geoutils::{Location, Distance};
//!
//! let berlin = Location::new(52.518611, 13.408056);
//! let moscow = Location::new(55.751667, 37.617778);
//! let is_in_radius = berlin.is_in_circle(&moscow, Distance::from_meters(2000.0)).unwrap();
//!
//! println!("Is Berlin in 2000m of Moscow? {}", is_in_radius);
//! ```
//!

#![deny(missing_docs)]
mod formula;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use formula::Distance;

/// Location defines a point using it's latitude and longitude.
#[derive(Debug, PartialEq, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
    /// This method returns Err if the formula fails to converge within 100 iterations.
    pub fn distance_to(&self, to: &Location) -> Result<Distance, String> {
        match formula::vincenty_inverse(self, to, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance),
            Err(e) => Err(e),
        }
    }

    /// Find the distance from itself to another point using Haversine formula.
    /// This is usually computationally less intensive than [distance_to](struct.Location.html#method.distance_to) but
    /// is generally not as accurate.
    pub fn haversine_distance_to(&self, to: &Location) -> Distance {
        formula::haversine_distance_to(self, to)
    }

    /// Check if the point is within a fixed radius of another point.
    pub fn is_in_circle(&self, center: &Location, radius: Distance) -> Result<bool, String> {
        match formula::vincenty_inverse(self, center, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance.meters() < radius.meters()),
            Err(e) => Err(e),
        }
    }

    /// Find the center of given locations.
    pub fn center(coords: &[&Location]) -> Location {
        formula::center_of_coords(coords)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_distance() {
        let l1 = Location::new(27.740068, 85.337576);
        let l2 = Location::new(27.740286, 85.337059);

        match l1.distance_to(&l2) {
            Ok(distance) => {
                assert_eq!(distance.meters(), 56.409);
            }
            Err(e) => panic!("Failed: {:?}", e),
        }
    }

    #[test]
    fn test_get_distance_haversine() {
        let l1 = Location::new(27.740068, 85.337576);
        let l2 = Location::new(27.740286, 85.337059);

        let distance = l1.haversine_distance_to(&l2);
        assert_eq!(distance.meters(), 56.36);
    }

    #[test]
    fn test_get_center() {
        let l1 = Location::new(52.518611, 13.408056);
        let l2 = Location::new(55.751667, 37.617778);

        let l = Location::center(&[&l1, &l2]);

        assert_eq!(l.0, 54.743683);
        assert_eq!(l.1, 25.033239);
    }
}
