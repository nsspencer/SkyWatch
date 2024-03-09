# Astro Access

Using astropy coordinates, this library allows you to compute access from one object in the galaxy to another.

## Constraints

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