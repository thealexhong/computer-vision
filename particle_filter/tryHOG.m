% Alex Hong, 2014
% Read the README file to see how to use this script file
% Try out the HOG detectors on a moused in patch and orientation
run '../../matlab/startup.m'

if ~exist('seqRoot')
  % YOU NEED TO CHANGE THIS
  seqRoot = '../pigeon/';
end

% YOU SHOULD RESET THESE TO TRY DIFFERENT TARGETS
pigMode = 2   % Use 1 or 2  (for contrasty pigeon or bland pigeon, resp)
marking = 1   % Use 0 or 1 (for back marking, or head, resp)

if marking == 0
  resDir = './back/';
elseif marking == 1
  resDir = './head/';
end

load(sprintf('%shog%dmark%d', resDir, pigMode, marking));
histModel = hogParams.histModel;
wCell = hogParams.wCell;

frameStart = 30;
frameEnd = 180;
frameNum = (frameStart-1)+ceil(rand(1,1)*(frameEnd - frameStart + 1));

im = pgmRead([seqRoot 'pigClone' num2strPad(frameNum, 4) '.pgm']);

figure(1); clf; showIm(im);
title(['Frame ' sprintf('%d', frameNum)]);

fprintf(2, 'Input bounding box (top left, bot right)...\n');
cb = round(ginput(2));
imCrop = im(cb(1,2):cb(2,2), cb(1,1):cb(2,1));

figure(1); clf; showIm(imCrop);
title('Test image');

fprintf(2,'Mouse in two points, p(1,:), p(2,:) such that their difference\n');
fprintf(2,'p(2,:)-p(1,:) is a direction vector to be rotated to be\n');
fprintf(2, 'horizontal to the right (i.e. ---->)\n');
cb = ginput(2);
dr = cb(2,:)-cb(1,:);
dr = dr'/norm(dr);
A = eye(3);
A(1:2, 1:2) = [dr [-dr(2); dr(1)]];
mn = (1+[size(imCrop,2); size(imCrop,1)])/2.0;
A(1:2, 3) = mn - A(1:2, 1:2)*mn;
% Do the warp
[imB, B] = imWarpAffine(imCrop, A, 1);
imB(isnan(imB)) = 0;

% Crop the result.
im1 = imB(1:size(imCrop,1), 1:size(imCrop,2));

figure(2); clf; showIm(im1);
title('Rotated test image');




[x,y] = meshgrid(1:2:size(im1,2), 1:2:size(im1,1));
pts = [x(:) y(:)]';

nPts = size(pts,2);
hist = getAlignedHOG(pts, im1, hogParams);
hist = reshape(hist, [hogParams.nyCell, hogParams.nxCell, ...
                    hogParams.nAmp * hogParams.nTheta, nPts]);

sumX2 = zeros(nPts,1);
cellResp = zeros(hogParams.nyCell, hogParams.nxCell, nPts);
for ix = 1:hogParams.nxCell
  for iy = 1:hogParams.nyCell
    tmp = sum(hist(iy,ix,:,:).* ...
              repmat(histModel(iy, ix,:), [1,1,1,nPts]),3)./ ...
          sqrt((sum(hist(iy,ix,:,:).^2,3)+1)*...
               (sum(histModel(iy,ix,:).^2)+1));
    cellResp(iy,ix,:) = tmp;
    sumX2 = sumX2 + reshape(wCell(iy, ix, :) * tmp, nPts, 1);
  end
end

resp = sumX2/sum(wCell(:));

figure(10); clf; showIm(reshape(resp, size(x)));
title('Histogram similarity image');

lgOdds = logOdds(resp, hogParams.lambda, hogParams.mu, hogParams.scl);
figure(11); clf; showIm(reshape(lgOdds, size(x)));
title('Log P-on / P-off image');

figure(12); clf; 
showIm( exp(reshape(lgOdds, size(x))));
title('P-on / P-off image');




