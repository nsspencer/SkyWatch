# SkyPath

SkyPath is an extension of Astropy's SkyCoord class which is useful for representing a body's coordinate states over time.

## Main features:

1) Fast position and velocity interpolations at any time(s) and coordinate frame.
2) Access functions for defining intervals of time when one body can "see" another.
3) Look angle (Azimuth, Elevation, and Range) calculations from one SkPath's custom defined body frame (with attitude) to another.


## Use Cases:

* Modeling position and velocity components of physical objects in any coordinate frame over time
* Satellite communication modeling
* Plantery movement simulations
* Aircraft tracking
* Celestial body observations


## Access

### Constraints

Access is an amorphous term, but this library allows you to define what it means to you through the use of *constraints*.

Constraints act on properties of the observer, target, and time to determine whether or not access is successful.

Typically, constraints represent physical properties such as position and velocity.

The following constraints are currently supported:

#### Line Of Sight

Determined by calculating when a body blocks visibility between two coordinates.

#### Azimuth, Elevation, and Range

When the azimuth, elevation, and range from an observer to a target are within a given min/max value.

#### Temporal

A time window which determines when access is successful.


## Line Of Sight (AER)

### Body Frames
* East, North, Up
* Nadir facing with velocity constraint
