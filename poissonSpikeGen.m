function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
dt=1e-3;
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;