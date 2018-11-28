use std::f64::consts::PI;

/// Location defines a point using it's latitude and longitude.
pub struct Location(f64, f64);

impl Location {
    pub fn new(lat: f64, lon: f64) -> Self {
        Location(lat, lon)
    }

    /// Find the distance from itself to another point.
    pub fn distance_to(self, to: Location) -> Result<f64, String> {
        match vincenty_inverse(self, to, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance),
            Err(e) => Err(e),
        }
    }

    /// Check if the point is within a fixed radius of another point.
    pub fn is_in_circle(self, center: Location, radius: f64) -> Result<bool, String> {
        match vincenty_inverse(self, center, 0.00001, 0.0) {
            Ok(res) => Ok(res.distance < radius),
            Err(e) => Err(e),
        }
    }

    // Find the center of given locations.
    pub fn center(coords: Vec<Location>) -> Location {
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

        Location(lat.to_degrees(), lon.to_degrees())
    }
}

/// DistanceResult is a private struct to hold results from Vincenty's inverse formula.
struct DistanceResult {
    distance: f64,
    _initial_bearing: Option<f64>,
    _final_bearing: Option<f64>,
    _iterations: i32,
}

/// Compute the Vincenty Inverse formula on two points 'start' and 'end'.
fn vincenty_inverse(
    start: Location,
    end: Location,
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
    let (mut sin_lambda, mut cos_lambda) = (0.0, 0.0);
    loop {
        sin_lambda = lambda.sin();
        cos_lambda = lambda.cos();
        let sin_sq_sigma = (cos_u2 * sin_lambda) * (cos_u2 * sin_lambda)
            + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda)
                * (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda);

        // Points coincide
        if sin_sq_sigma == 0.0 {
            break;
        }

        sin_sigma = sin_sq_sigma.sqrt();
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
        sigma = sin_sigma.atan2(cos_sigma);
        let sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma;
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
    let A = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let B = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));

    let delta_sigma = B
        * sin_sigma
        * (cos2_sigma_m
            + B / 4.0
                * (cos_sigma * (-1.0 + 2.0 * cos2_sigma_m * cos2_sigma_m)
                    - B / 6.0
                        * cos2_sigma_m
                        * (-3.0 + 4.0 * sin_sigma * sin_sigma)
                        * (-3.0 + 4.0 * cos2_sigma_m * cos2_sigma_m)));
    let s = b * A * (sigma - delta_sigma);

    let mut alpha_1 = cos_u2 * sin_lambda.atan2(cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda);
    let mut alpha_2 = cos_u1 * sin_lambda.atan2(-sin_u1 * cos_u2 + cos_u1 * sin_u2 * cos_lambda);

    alpha_1 = (alpha_1 + 2.0 * PI) % (2.0 * PI);
    alpha_2 = (alpha_2 + 2.0 * PI) % (2.0 * PI);

    Ok(DistanceResult {
        distance: (s * 1000.0).round() / 1000.0,
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

        match l1.distance_to(l2) {
            Ok(distance) => {
                assert_eq!(distance, 56.409);
            }
            Err(e) => panic!("Failed: {:?}", e),
        }
    }

    #[test]
    fn test_get_center() {
        let l1 = Location::new(52.518611, 13.408056);
        let l2 = Location::new(55.751667, 37.617778);

        let l = Location::center(vec![l1, l2]);

        assert_eq!(l.0, 54.743683);
        assert_eq!(l.1, 25.033239);
    }
}
