# SkyPath

SkyPath is an extension of Astropy's SkyCoord class which is useful for representing a body's coordinate states over time.

## Main features:

1) Fast position and velocity interpolations in any astropy coordinate frame.
2) Access functions for defining intervals of time when multiple user defined (or builtin) constraints are satisfied.
3) Look angle (Azimuth, Elevation, and Range) calculations from one coordinate to another.


## Use Cases:

* Modeling position and velocity components of physical objects in any coordinate frame over time
* Satellite communication modeling
* Plantery movement simulations
* Aircraft tracking
* Celestial body observations


## Examples:

See tests\test_smoke.py and tests\test_look_angles.py
