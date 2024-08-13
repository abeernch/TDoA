CEP50 Computation and plotting.

Main script: 
CEP50_220523

Main Functions:
1. computeCEP50
2. computeCRLB

Auxiliary Functions:
1. jacobian
2. ensureInvertible
3. constants (class)
4. parseReferenceSensor
5. resampleCovMtx
6. excludeFromLegend

CONFIGURABLE PARAMETERS:

1. nSensors: No. of sensors in the network 
2. timingError: To create the error covariance matrix which is used to compute the CRLB
3. grd_size: Resolution of grid representing potential positions of source
4. x_range: grid limits/region over which to analyze network performance along x-axis
5. y_range: grid limits/region over which to analyze network performance along y-axis