function newSamples = resampleRandom(sampleDist)

  numSamps = size(sampleDist.samples,1);

  % Build Cumulative Distribution
  cdf = [0; cumsum(sampleDist.weights)];
  
  x = rand(numSamps,1);
  % Sample from the source weighted (discrete) distribution
  for n = 1:numSamps
    m = find(cdf >= x(n));
    cdfindex = m(1)-1;
    newSamples.samples(n,:) = sampleDist.samples(cdfindex,:);
  end
  newSamples.weights = ones(numSamps,1) / numSamps;





