use super::Location;
use std::f64::consts::PI;
use std::fmt;

/// Distance represents a physical distance in a certain unit.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Distance(f64);

impl fmt::Display for Distance {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} meters", self.meters())
    }
}

impl Distance {
    /// Create a distance in meters.
    pub fn from_meters<P: Into<f64>>(m: P) -> Self {
        return Distance(m.into());
    }

    /// Return the distance in meters.
    pub fn meters(&self) -> f64 {
        return self.0;
    }
}

/// DistanceResult is a struct to hold results from Vincenty's inverse formula.
pub struct DistanceResult {
    pub distance: Distance,
    _initial_bearing: Option<f64>,
    _final_bearing: Option<f64>,
    _iterations: i32,
}

/// Compute the Vincenty Inverse formula on two points 'start' and 'end'.
pub fn vincenty_inverse(
    start: &Location,
    end: &Location,
    _accuracy: f64,
    _precision: f64,
) -> Result<DistanceResult, String> {
    // WGS-84 geocentric datum parameters
    let a: f64 = 6378137.0; // Semi-major axis
    let b: f64 = 6356752.314245; // Semi-minor axis
    let f: f64 = 1.0 / 298.257223563; // Inverse-flattening

    // Start and end points in Radians
    let p1 = (start.0.to_radians(), start.1.to_radians());
    let p2 = (end.0.to_radians(), end.1.to_radians());

    // Difference in longitudes
    let l = p2.1 - p1.1;

    // u = 'reduced latitude'
    let (tan_u1, tan_u2) = ((1.0 - f) * p1.0.tan(), (1.0 - f) * p2.0.tan());
    let (cos_u1, cos_u2) = (
        1.0 / (1.0 + tan_u1 * tan_u1).sqrt(),
        1.0 / (1.0 + tan_u2 * tan_u2).sqrt(),
    );
    let (sin_u1, sin_u2) = (tan_u1 * cos_u1, tan_u2 * cos_u2);

    // First approximation
    let mut lambda = l;
    let mut iter_limit = 100;
    let mut cos_sq_alpha = 0.0;
    let (mut sin_sigma, mut cos_sigma, mut cos2_sigma_m, mut sigma) = (0.0, 0.0, 0.0, 0.0);
    let (mut _sin_lambda, mut _cos_lambda) = (0.0, 0.0);
    loop {
        _sin_lambda = lambda.sin();
        _cos_lambda = lambda.cos();
        let sin_sq_sigma = (cos_u2 * _sin_lambda) * (cos_u2 * _sin_lambda)
            + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * _cos_lambda)
                * (cos_u1 * sin_u2 - sin_u1 * cos_u2 * _cos_lambda);

        // Points coincide
        if sin_sq_sigma == 0.0 {
            break;
        }

        sin_sigma = sin_sq_sigma.sqrt();
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * _cos_lambda;
        sigma = sin_sigma.atan2(cos_sigma);
        let sin_alpha = cos_u1 * cos_u2 * _sin_lambda / sin_sigma;
        cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
        cos2_sigma_m = if cos_sq_alpha != 0.0 {
            cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_sq_alpha
        } else {
            0.0
        };
        let c = f / 16.0 * cos_sq_alpha * (4.0 + f * (4.0 - 3.0 * cos_sq_alpha));
        let lambda_prime = lambda;
        lambda = l
            + (1.0 - c)
                * f
                * sin_alpha
                * (sigma
                    + c * sin_sigma
                        * (cos2_sigma_m
                            + c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m * cos2_sigma_m)));

        iter_limit -= 1;
        if (lambda - lambda_prime).abs() > 1e-12 && iter_limit > 0 {
            continue;
        }

        break;
    }

    if iter_limit <= 0 {
        return Err("formula failed to converge".to_string());
    }

    let u_sq = cos_sq_alpha * (a * a - b * b) / (b * b);
    let cap_a = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let cap_b = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));

    let delta_sigma = cap_b
        * sin_sigma
        * (cos2_sigma_m
            + cap_b / 4.0
                * (cos_sigma * (-1.0 + 2.0 * cos2_sigma_m * cos2_sigma_m)
                    - cap_b / 6.0
                        * cos2_sigma_m
                        * (-3.0 + 4.0 * sin_sigma * sin_sigma)
                        * (-3.0 + 4.0 * cos2_sigma_m * cos2_sigma_m)));
    let s = b * cap_a * (sigma - delta_sigma);

    let mut alpha_1 = cos_u2 * _sin_lambda.atan2(cos_u1 * sin_u2 - sin_u1 * cos_u2 * _cos_lambda);
    let mut alpha_2 = cos_u1 * _sin_lambda.atan2(-sin_u1 * cos_u2 + cos_u1 * sin_u2 * _cos_lambda);

    alpha_1 = (alpha_1 + 2.0 * PI) % (2.0 * PI);
    alpha_2 = (alpha_2 + 2.0 * PI) % (2.0 * PI);

    Ok(DistanceResult {
        distance: Distance::from_meters((s * 1000.0).round() / 1000.0),
        _initial_bearing: if s == 0.0 {
            None
        } else {
            Some(alpha_1.to_degrees())
        },
        _final_bearing: if s == 0.0 {
            None
        } else {
            Some(alpha_2.to_degrees())
        },
        _iterations: iter_limit,
    })
}

/// Get the center of a list of coordinates.
pub fn center_of_coords(coords: &[&Location]) -> Location {
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);

    for loc in coords.iter() {
        let lat = loc.0.to_radians();
        let lon = loc.1.to_radians();

        x += lat.cos() * lon.cos();
        y += lat.cos() * lon.sin();
        z += lat.sin();
    }

    let number_of_locations = coords.len() as f64;
    x /= number_of_locations;
    y /= number_of_locations;
    z /= number_of_locations;

    let hyp = (x * x + y * y).sqrt();
    let lon = y.atan2(x);
    let lat = z.atan2(hyp);

    // Convert to degrees and round to 6 digits
    let lon = lon.to_degrees();
    let lon = (lon * 10.0_f64.powi(6)).round() / 10.0_f64.powi(6);

    let lat = lat.to_degrees();
    let lat = (lat * 10.0_f64.powi(6)).round() / 10.0_f64.powi(6);

    Location(lat, lon)
}

/// Implementation of Haversine distance between two points.
pub fn haversine_distance_to(start: &Location, end: &Location) -> Distance {
    let haversine_fn = |theta: f64| (1.0 - theta.cos()) / 2.0;

    let phi1 = start.latitude().to_radians();
    let phi2 = end.latitude().to_radians();
    let lambda1 = start.longitude().to_radians();
    let lambda2 = end.longitude().to_radians();

    let hav_delta_phi = haversine_fn(phi2 - phi1);
    let hav_delta_lambda = phi1.cos() * phi2.cos() * haversine_fn(lambda2 - lambda1);
    let total_delta = hav_delta_phi + hav_delta_lambda;

    let dist = (2.0 * 6371e3 * total_delta.sqrt().asin() * 1000.0).round() / 1000.0;
    Distance::from_meters(dist)
}
