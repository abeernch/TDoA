Created: 220928
Last edited: 220928

1. The script is compatible with any no of sensors and emitters
2. 2D TDoA solutions and isochrones
3. Compatible with time-varying (updating) emitters
4. Multiple solutions at the output. For example, in 2D, the solution may be obtained by the use of 03 sensors, anything above that is redundant. However, adding additional
sensors into the network allow for a more robust precision by ensuring maximum non-collinearity between emitter(s) and sensors. 
n+1 sensors are needed for solutions in n-dimensional TDoA scenarios. if no of sensors is > n+1, different baselines for a centralised sensor pairing can be formed to provide with multiple solutions. These multiple solutions for different baselines can be then used to eliminate non-unique and incorrect solutions by comparison. This features is 
a part of this version of the TDoA Algorithm