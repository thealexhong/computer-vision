% Alex Hong, 2014
function newSamples = updateWeights(sampDist, im, likelihoodParams)
  numSamps = size(sampDist.samples,1);
  
  % For each sample, first obtain data log likelihood
  loglike = pigLogLike(sampDist.samples', im, ...
                               likelihoodParams);
  
  % c = 1 / sum (exp(loglike));
  newSamples.weights = exp(loglike) / sum (exp(loglike));
  newSamples.samples = sampDist.samples;
  N = 1 / sum(newSamples.weights.^2);
  fprintf('Effective Number of Particles: %f\n', N);
  
  return;
