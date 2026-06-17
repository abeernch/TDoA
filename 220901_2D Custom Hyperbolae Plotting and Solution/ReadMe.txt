This folder contains early MATLAB work on 2D hyperbola construction and
coordinate transformation for TDoA-style isochrone plots.

Main functions:

1. `draw_2Disochrone` - generates the generic hyperbola and computes the
   parameters needed for coordinate transformation.
2. `coord_2Dtfm` - rotates, offsets, and recenters the generic hyperbola into
   the working reference frame.
3. `toa_generator` - generates TOA measurements for a sensor network and source
   configuration.

The scripts in this folder are experimental and represent one stage in the
development of the larger TDoA analysis workflow.
