[![Build Status](https://travis-ci.org/srishanbhattarai/geoutils.svg?branch=master)](https://travis-ci.org/srishanbhattarai/geoutils)
[![Documentation](https://docs.rs/geoutils/badge.svg)](https://docs.rs/geoutils/)

# geoutils

Geoutils is a evolving crate to provide several geological computations and utilities.
Most computations are based off methods on the [Location](struct.Location.html) struct.

Find the full API reference at [docs.rs](http://docs.rs/geoutils).

## Examples

* Get the distance between two points using [Vincenty's Inverse Formula](https://en.wikipedia.org/wiki/Vincenty%27s_formulae).
```rust
extern crate geoutils;

use geoutils::Location;

let berlin = Location::new(52.518611, 13.408056);
let moscow = Location::new(55.751667, 37.617778);
let distance = berlin.distance_to(&moscow).unwrap();

println!("Distance = {}", distance);
```

* Get the distance between two points using the [Haversine Formula](https://en.wikipedia.org/wiki/Haversine_formula).
```rust
extern crate geoutils;

use geoutils::Location;

let berlin = Location::new(52.518611, 13.408056);
let moscow = Location::new(55.751667, 37.617778);
let distance = berlin.haversine_distance_to(&moscow);

println!("Distance = {}", distance);
```

* Get the center of a list of coordinates.
```rust
extern crate geoutils;

use geoutils::Location;

let berlin = Location::new(52.518611, 13.408056);
let moscow = Location::new(55.751667, 37.617778);
let center = Location::center(vec![&berlin, &moscow]);

println!("Center {}, {}", center.latitude(), center.longitude());
```

* Check if a point falls in a certain radius of another point.
```rust
extern crate geoutils;

use geoutils::Location;

let berlin = Location::new(52.518611, 13.408056);
let moscow = Location::new(55.751667, 37.617778);
let is_in_radius = berlin.is_in_circle(&moscow, 2000);

println!("Is Berlin in 2000m of Moscow? {}", is_in_radius);
```


# License
Apache-2.0
