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

    /// Check if the point is the vertex of polygon or line given as vertices
    pub fn is_in_vertex(&self, vertices: &[&Location]) -> bool {
        for vertex in vertices {
            if vertex.latitude() == self.latitude() && vertex.longitude() == self.longitude() {
                return true;
            }
        }
        false
    }

    /// Check if the point is the edge of polygon or line given as vertices
    pub fn is_in_edge(&self, vertices: &[&Location]) -> bool {
        for (index, vertex) in vertices.iter().enumerate() {
            let mut next_index = index + 1;
            if index == vertices.len() - 1 {
                next_index = 0;
            }
            let next_vertex = vertices[next_index];
            if self.longitude() <= vertex.longitude().max(next_vertex.longitude())
                && self.longitude() >= vertex.longitude().min(next_vertex.longitude())
                && self.latitude() <= vertex.latitude().max(next_vertex.latitude())
                && self.latitude() >= vertex.latitude().min(next_vertex.latitude())
            {
                return true;
            }
        }
        false
    }
    /// Check if the point lies in polygon or line given as vertices
    pub fn is_contained_in(&self, vertices: &[&Location]) -> bool {
        let mut count: u32 = 0;
        let line1 = [self, &Location(f64::MAX, self.longitude())];
        for (index, vertex) in vertices.iter().enumerate() {
            let next_index = (index + 1) % vertices.len();
            let next_vertex = vertices[next_index];
            let line2 = [vertex, next_vertex];
            if formula::do_intersect(&line1, &line2) {
                count += 1;
            }
        }
        count % 2 == 1
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
    #[test]
    fn test_is_contained_in() {
        let colorado = [
            &Location(40.996549, -109.001960),
            &Location(40.983031, -102.125103),
            &Location(37.061297, -102.107194),
            &Location(37.075586, -109.001960),
        ];
        let colorado_corner = Location(37.075586, -109.001960);
        let denver = Location(39.725186, -104.984211);
        let chicago = Location(41.839142, -87.702715);
        assert_eq!(denver.is_contained_in(&colorado), true);
        assert_eq!(chicago.is_contained_in(&colorado), false);
        assert_eq!(colorado_corner.is_contained_in(&colorado), true);
    }
}
