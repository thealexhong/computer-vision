% Alex Hong, 2014

function newSamples = propagateSamples(sampDist, dynamicsParams)
  sigma = diag(dynamicsParams.dynCovar);
  N = size(sampDist.samples, 1);
  newSamples.samples = normrnd(sampDist.samples, ...
                               repmat(sigma', N, 1));
  newSamples.weights = sampDist.weights;
end

