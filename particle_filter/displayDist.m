function im = displayDist(sampDist, s)

  im = zeros(s);
  numSamps = size(sampDist.samples,1);
  samps = sampDist.samples;
  w = sampDist.weights;
  
  for n = 1:numSamps
      x = round(samps(n,1:2));
      x = min( [x; s(2:-1:1)] );
      x = max( [x; 1 1] );
      im(x(2), x(1)) = im(x(2), x(1)) + w(n);
  end
  
  sigmaBlur = 1.0;
  sigmaSqr = sigmaBlur*sigmaBlur;
  gFiltSize = 2 * round(3.0 * sigmaBlur) + 1;
  x = [1:gFiltSize] - round((gFiltSize+1)/2);
  gFilt = exp(- x .* x / (2.0*sigmaSqr));
  gFilt = gFilt/ sum(gFilt(:));

  im = rconv2sep(im, gFilt, gFilt);
  return;
