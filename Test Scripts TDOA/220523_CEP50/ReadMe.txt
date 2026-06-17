CEP50 computation and plotting.

Main script:
`CEP50_220523`

Main functions:
1. `computeCEP50`
2. `computeCRLB`

Auxiliary functions:
1. `jacobian`
2. `ensureInvertible`
3. `constants` (class)
4. `parseReferenceSensor`
5. `resampleCovMtx`
6. `excludeFromLegend`

Configurable parameters:
1. `nSensors` - number of sensors in the network.
2. `timingError` - error covariance used for CRLB computation.
3. `grd_size` - grid resolution for source-position analysis.
4. `x_range` - x-axis extent for the evaluation grid.
5. `y_range` - y-axis extent for the evaluation grid.
