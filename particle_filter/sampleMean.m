function m = sampleMean(sampDist)

  w = sampDist.weights;
  s = sampDist.samples;

  m = sum(s.*repmat(w,1,size(s,2)));
