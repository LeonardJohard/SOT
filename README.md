# SOT
Matlab code for self-organizing trajectories

Takes a set of n-dimensional trajectories and produces a  piecwise linear 'mean' trajectory. Calculates intrastroke variance around each linesegment.

Standard SOM can be viewed as a 2-step process:
1. Assign each point to the closest node
2. Adapt position of closest node and its neighbors towards the mean of the assigned points

SOT works on a set of points as follows:
1. Calculate for each trajectory (a set of consecutive points) the closest DTW mapping to the saved SOM trajectory (a set of consecutive nodes).
2. Assign each point to the corresponding node of the DTW mapping.
3. As 2 above

SOT is a restriction of the regular SOM point to keep the consecutive order of the point-to-node mapping.


Usage: Add() all trajectories, then run adapt() and adaptDTW() until
convergence. Get variance from segvar estimate or run fullmean for
complete analysis.


See example.m and example2.m for rowing trajectory examples

Note that this implementation is based on an recovered earlier 2013 version of the code and might differ in details from the implementation in the paper. I have tried to recreate the code of the paper to the best of my ability.

The implementation uses the Dynamic Time Warp implementation by Dan Ellis (http://www.ee.columbia.edu/~dpwe/resources/matlab/dtw/), which needs to be compiled.
