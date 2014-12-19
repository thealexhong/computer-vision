% Modified by Alex Hong, 2014

%% FILE: browseResults.m
%% Simple demonstration of a particle filter for visual tracking
%%
%% Authors: ADJ and DJF, Feb. 2006, Modified Nov. 2007

clear all
%close all
run '../../matlab/startup.m'
question3 = true;

if (question3)
    % YOU NEED to change these for Question 3
    cropDisplay = false;  % Set to false for question 3
    outPrefix = 'pig';   % prefix of mat file for tracking results
                     % Eg. partFilt writes to results/pig#mark#frm#.mat
                     % so the prefix is 'pig'
else
    % YOU NEED to change these for Question 3
    cropDisplay = true;  % Set to false for question 3
    outPrefix = 'pig';   % prefix of mat file for tracking results
                     % Eg. partFilt writes to results/pig#mark#frm#.mat
                     % so the prefix is 'pig'
end
    cropScl = 2;         % Blow up the result image Figure 2 by
                     % this factor for display and writing to jpeg.

for pigMode = 3;  % 1 or 2 or 3
  
  marking = 1;  % 0 or 1



  %  Read sequence and display
  seqParams.frameRoot = '../pigeon/';
  seqParams.framePrefix = 'pigClone';
  seqParams.padNum = 4;
  seqParams.frameExt = 'ppm';
  seqParams.frameStart = 30;
  seqParams.frameSkip = 2;
  seqParams.frameEnd = 180;
  resDir = './results/';

  % enlargement factor for display
  seqParams.figScl = 1;
  if marking == 0
    mn2MarkDist = 0;
  elseif marking == 1
    mn2MarkDist = 63;
  end
  sqrCropRad = 75;

  load 'clonePigCoarseRes';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % parameters for state space definition

  dim = 3;

  numsamples = 1000;

  % Loop over the image filenames in fn and go nuts
  frameNum0 = seqParams.frameStart;
  for frameNum = frameNum0:seqParams.frameSkip:seqParams.frameEnd
    % Read and display the frameNum^th image
    im = pgmRead([seqParams.frameRoot seqParams.framePrefix ...
                  num2strPad(frameNum, 4) '.pgm']);

    figure(1); clf; showIm(im);
    title(sprintf('Frame %d', frameNum));
    pause(0.1);
    
    if cropDisplay
      % Crop rough neighbourhood of target region.
      bodyMn = coarseRes{frameNum}.mn(pigMode,:);
      U = coarseRes{frameNum}.U(:,:,pigMode);
      
      xc0 = round(bodyMn + U(:,1)' * mn2MarkDist);
      xRange = xc0(1) + ((-sqrCropRad):sqrCropRad);
      yRange = xc0(2) + ((-sqrCropRad):sqrCropRad);
      xRange = xRange(xRange >=1 & xRange <= size(im,2));
      yRange = yRange(yRange >=1 & yRange <= size(im,1));

      imCrop = im(yRange, xRange);
    
      figure(1); hold on;
      plot([xRange(1) xRange(1) xRange(end) xRange(end) xRange(1)],...
           [yRange(1) yRange(end) yRange(end) yRange(1) yRange(1)], 'g');
      spaceGrid = 25;
    else
      imCrop = im;
      xRange = [1, size(im,2)];
      yRange = [1, size(im,1)];
      spaceGrid = round(max(xRange(2), yRange(2))/15);
    end

    %% Get sampDist results
    load(sprintf('%s%s%dmark%dfrm%s.mat', resDir, outPrefix, ...
                 pigMode, marking, ...
                 num2strPad(frameNum, seqParams.padNum)));
    nSamp = size(sampDist.samples, 1);
    sampDist.samples = sampDist.samples - ...
        repmat([xRange(1), yRange(1), 1]-1, nSamp,1);

    %% Display image and results
    
    im = zeros(2*size(imCrop,1), size(imCrop,2));
    im(1:size(imCrop,1),:) = imCrop(:,:)/255;
    
    lm = -log(displayDist(sampDist,size(imCrop(:,:)))+1/numsamples);
    lm = lm - min(lm(:));
    lm = lm/max(lm(:));
    
    im((size(imCrop,1)+1):end,:) = lm;
    
    % Draw grid lines on image
    im(1:spaceGrid:size(imCrop,1), :) = 0;
    im(size(imCrop,1)+(1:spaceGrid:size(imCrop,1)), :) = 0;
    im(:, 1:spaceGrid:size(imCrop,2)) = 0;
    
    % Draw orientation posterior
    thetas = sampDist.samples(:,3);
    idx = thetas >= 2*pi;
    while any(idx)
      thetas(idx) = thetas(idx) - 2*pi;
      idx = thetas > 2*pi;
    end
    idx = thetas < 0;
    while any(idx)
      thetas(idx) = thetas(idx) + 2*pi;
      idx = thetas < 0;
    end
    
    nHist = 101;
    thetaHist = zeros(nHist,1);
    thetaBin = (0:(nHist-1))*(2*pi/nHist);
    for k=1:nHist
      dist = thetas - thetaBin(k);
      idx = dist >= pi;
      while any(idx)
        dist(idx) = dist(idx)- 2*pi;
        idx = dist >= pi;
      end
      idx = dist < -pi;
      while any(idx)
        dist(idx) = dist(idx) + 2*pi;
        idx = dist < -pi;
      end
      idx = -pi/nHist <= dist & dist < pi/nHist;
      if any(idx)
        thetaHist(k) = sum(sampDist.weights(idx));
      end
    end
    figure(3); clf;
    plot(thetaBin, thetaHist/sum(thetaHist), 'LineWidth',2);
    ax = axis;
    axis([0 2*pi 0, ax(4)]);
    set(gca, 'FontSize', 14);
    xlabel('Particle angle: \theta');
    ylabel('Relative frequency');
    title('Orientation distribution');
    
    
    meanPosterior = sampleMean(sampDist);
    mxy = meanPosterior(1:2);
    r = mean([cos(thetas) sin(thetas)], 1);
    theta = atan2(r(2), r(1));
    tang = [cos(theta), sin(theta)];
    tang = tang; nrml = [tang(2), -tang(1)];
    poly = [mxy + 3*nrml;
            mxy + 10*tang;
            mxy - 3*nrml ];
    %  template = sgnDist2Poly(1:size(im,2), 1:size(im,1), poly);  
    %  im(template<0) = 1;
    % Add black border on top and left sides.
    im0 = zeros(1+size(im));
    im0(2:end, 2:end) = im;
    figure(2); clf;
    showIm(im0);  
    resizeImageFig(figure(2), size(im), cropScl); 
    hold on;
    h = plot(1+poly(:,1),1+poly(:,2), 'r', 'LineWidth', 2);
    h = plot(1+mxy(1),1+mxy(2), 'or', 'MarkerSize', 5, ...
             'MarkerFaceColor', 'r');
    
    % Remember to set page setup on figure window to match screen size.
    %% Get sampDist results
    fOut = sprintf('%s%s%dmark%dfrm%s', resDir, outPrefix, pigMode, ...
                   marking, num2strPad(frameNum, seqParams.padNum));
    print(figure(2), '-djpeg95', fOut);
    
    pause(0.1);  %% ensure plots are drawn
  end
end
