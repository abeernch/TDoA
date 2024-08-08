1. This folder contains work on the custom algorithm created for the formation of hyperbolae as isochrones for the TDoA position estimation algorithm

FUNCTIONS:

1. draw_2Disochrone: (CREATES GENERIC HYPERBOLA AND COMPUTES PARAMETERS FOR COORD TFM)
This function creates hyperbolae in the opening up/down configuration. The inputs to the function are the sensor network geometry parameters, 4-quad angle of the slave 
sensor w.r.t the reference sensor and the RDOA with respect to the reference sensor . This function is called iteratively for n tdoa pairs to plot n hyperbolae.

2. coord_2Dtfm: (PERFORMS COORDINATE TRANSFORMATION)
- This is a nested function that is called within the isochrone_func. 
- Once a generic hyperbola is created, the necassary coordinate transformation from body to fixed frame
reference is completed here  w.r.t the reference sensor. 
- The transformation angle, adjusted offset center of the hyperbola and RDOA are input. 
- The Rotation matrix is calculated and matrix multiplied with the coordinates of the hyperbolae
- Based on the sign of RDOA, the calculated offset is either added (-ive rdoa) or subtracted (+ive rdoa) from the rotated points 

- Once the transformation is complete, the correct leaf of the hyperbola is to be selected. This is done again by keeping sign of the RDOA as a check. Based on the 
sign of the RDOA one of two leaves are selected.
- Lastly, to re-center the hyperbola to the reference sensor the position coords of the reference sensor are added to the coords of the hyperbola

3. toa_generator: (GENERATES TOA MEASUREMENTS FOR GIVEN SENSOR NETWORK AND SOURCE PARAMETERS)
- The positions for the sensors and the emitter can be defined interactively within the function, or they can be manually input if the interactive plotting is 
disabled by setting the 'graph' input as false
- Other inputs are total time and PRI. These are used to create an array of TOAs based on a preselected PRI for all sensors for the entered total time
- Based on the interactively (or otherwise) placed platforms, the ranges are calculated and converted to TOA at each sensor. 
- The is adapated to cater for any no of sensors and emitters 
- At the output, the platform positions, toa array and no of platforms are acquired.