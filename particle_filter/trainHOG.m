% Alex Hong, 2014
% See the README file for how to use this script.
run '../../matlab/startup.m'

if ~exist('seqRoot')
  % CHANGE THIS
  seqRoot = '../pigeon/';
end

% RESET THESE :
pigMode = 2 % 1  (Textured pigeons), 2 (bland pigeon).
pOn = false  % true or false
marking = 1 % 0  for back, or 1 for head

% The rest does not need to be changed.
if (marking == 0)
  resDir = './back/';
  mn2MarkDist = 0; % Used to crop rough marking region for training pOn.
elseif (marking == 1)
  resDir = './head/';
  mn2MarkDist =  63;
end
sqrCropRad = 100;

% Coarse tracking results.
load 'clonePigCoarseRes'

load(sprintf('%shog%dmark%d', resDir, pigMode, marking));
histModel = hogParams.histModel;
wCell = hogParams.wCell;

nRespHist = 101;
respHist = zeros(nRespHist,1);
frameStart = 30;
frameEnd = 180;
while (true)  % Hit Enter without mousing in a point to break out
  
  frameNum = (frameStart-1)+ceil(rand(1,1)*(frameEnd - frameStart + 1));
  frameNum = 2*floor(frameNum/2);  % Use only even frames.
  im = pgmRead([seqRoot 'pigClone' num2strPad(frameNum, 4) '.pgm']);

  if pOn
    % Crop rough neighbourhood of mark region.
    bodyMn = coarseRes{frameNum}.mn(pigMode,:);
    U = coarseRes{frameNum}.U(:,:,pigMode);
    
    xc0 = round(bodyMn + U(:,1)' * mn2MarkDist);
    xRange = xc0(1) + ((-sqrCropRad):sqrCropRad);
    yRange = xc0(2) + ((-sqrCropRad):sqrCropRad);
    xRange = xRange(xRange >=1 & xRange <= size(im,2));
    yRange = yRange(yRange >=1 & yRange <= size(im,1));

    imCrop = im(yRange, xRange);
    im = imCrop;    
  end
  
  figure(1); clf; showIm(im);
  title('Original image');

  if pOn
    fprintf(2, 'Mouse in center of marking...\n');
    xc = round(ginput(1));
    if size(xc,2) == 0
      break;
    end
    [x, y] = meshgrid(xc(1)+(-3:3), xc(2)+(-3:3));
    imCrop = im;
    fprintf(2, 'Input one point in front of the marking...\n');
    cb = ginput(1);
    if size(cb,2) == 0
      break;
    end
    cb = [xc; cb];
  else
    fprintf(2, 'Input bounding box...\n');
    cb = round(ginput(2));
    if size(cb, 2) == 0
      break;
    end
    imCrop = im(cb(1,2):cb(2,2), cb(1,1):cb(2,1));
    
    figure(1); clf; showIm(imCrop);
    title('Cropped image');
    fprintf(2, 'Input two points indicating orientation ...\n');
    cb = ginput(2);
    if size(cb, 2) == 0
      break;
    end
  end
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
  title('Cropped and rotated image');
  if pOn
    nPts = prod(size(x));
    pts = inv(A) * [x(:), y(:), ones(nPts,1)]';
    pts = pts(1:2, :);
  else
    [x,y] = meshgrid(1:2:size(im1,2), 1:2:size(im1,1)); 
    pts = [x(:) y(:)]';
  end

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
  title('Histogram Similarity');

  val = ceil(resp * nRespHist);
  for k = 1:nRespHist;
    respHist(k) = respHist(k) + sum(val==k);
  end
  
  figure(100); clf; plot(respHist,'b');
  if pOn
    title('Histogram for P-on');
  else
    title('Histogram for P-off');
  end
  
  figure(101); clf; plot(log(respHist + 0.1),'b');
  if pOn
    title('Log-Histogram for P-on');
  else
    title('Log-Histogram for P-off');
  end
  
  fprintf(2, 'Press any key to continue...');
  pause;
  fprintf(2, '\n');
end

if true
  % Save resp histogram
  if pOn
    save(sprintf('%srespHistPig%dmark%dOn', resDir, pigMode, marking),...
       '-V6', 'respHist');
  else
    save(sprintf('%srespHistPig%dmark%dOff', resDir, pigMode, marking),...
       '-V6', 'respHist');
  end
end
  