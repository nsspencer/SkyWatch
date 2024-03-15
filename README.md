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

Access is an amorphous term, but this library allows you to define what it means to you through the use of *constraints*.

Constraints act on properties of the observer, target, and time to determine whether or not access is successful.

Typically, constraints represent physical properties such as position and velocity.

The following constraints are currently supported:

### Line Of Sight

Determined by calculating when a body blocks visibility between two coordinates.

### Azimuth, Elevation, and Range

When the azimuth, elevation, and range from an observer to a target are within a given min/max value.

### Temporal

A time window which determines when access is successful.


## Look Angles

Look angles are defined in Azimuth, Elevation, and Rotation in reference to a SkyCoord's body frame.

There are different strategies for calculating look angles from one SkyCoord to another.

## East North Up

This frame is useful for terresetrial based SkyPaths, where the east vector faces longitudinal east, north faces true north, and the up
vector the zenith projection directly out of the plane tangent to the terrestrial surface.


## Nadir Facing with Velocity Constraint

This frame is useful for platforms that have a velocity component and some attitude rotation about their body frame. Satellites and airplanes are likely
to be defined in this frame. The X vector is the direction of forward velocity, Z is nadir to the center of the coordinate frame, and Y is the cross
product of X and Z to complete the right handed system.

