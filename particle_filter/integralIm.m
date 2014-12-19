function intIm = integralIm(im)
% intIm = integralIm(im)

intIm = cumsum(im,2);
intIm = cumsum(intIm,1);

return;
