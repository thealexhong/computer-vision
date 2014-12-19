% Modified by Alex Hong, 2014

%% FILE: partFilt.m
% Simple demonstration of a particle filter for visual tracking
%
% Authors: ADJ and DJF, Feb. 2006  Revised Nov. 2007.

% clear
% close all

run '../../matlab/startup.m'
question3 = true;

if ~exist('seqRoot')
  % CHANGE THIS
  seqRoot = '../pigeon/';
end
if (question3)
    % For question 2 leave these flags alone
    safetyNet = false;      % Use coarse image to crop a rough image
                       % region around the target.
    cropRotatedIm = false;  % Ignore the corners when rotating the
                       % cropped image.
    sqrCropRad = 400;       % Manhattan radius of the cropped image.

    outPrefix = 'pig';   % Prefix of output mat file for tracking results
                     % Eg. partFilt writes to results/pig#mark#frm#.mat
                     % with the prefix 'pig'
else
    % For question 2 leave these flags alone
    safetyNet = true;      % Use coarse image to crop a rough image
                       % region around the target.
    cropRotatedIm = true;  % Ignore the corners when rotating the
                       % cropped image.
    sqrCropRad = 75;       % Manhattan radius of the cropped image.

    outPrefix = 'pig';   % Prefix of output mat file for tracking results
                     % Eg. partFilt writes to results/pig#mark#frm#.mat
                     % with the prefix 'pig'
end
                     
cropScl = 2;  % Blow up the cropped image by this factor for display.

% YOU NEED to reset pigMode and marking to track different targets
for pigMode = 3; % 1 denotes grey pigeon, 2 denotes red pigeon, 3 clone
                   
  marking = 1;  % 0 for back, 1 for head

  % Read sequence and display
  seqParams.frameRoot = seqRoot;
  seqParams.framePrefix = 'pigClone';
  seqParams.padNum = 4;
  seqParams.frameExt = 'pgm';
  seqParams.frameStart = 30;
  seqParams.frameSkip = 2;
  seqParams.frameEnd = 180;

  % enlargement factor for display
  seqParams.figScl = 1;

  % Load coarsely tracked results.
  load 'clonePigCoarseRes';

  % Initialization
  if marking == 0
    resDir = './back/';
    mn2HeadDist = 0;
    % Mean initial guess for first frame.
    meanPos = [529   145  -135/180*pi;...
               159   323 35/180*pi;...
               272   180 135/180*pi];
    if (pigMode == 2)
        meanPos = repmat([159, 323, 35 / 180*pi], 3, 1);
    end
           
           
  elseif marking == 1
    resDir = './head/';
    mn2HeadDist = 63;
    % Mean initial guess for first frame.
    meanPos = [470   116  -160/180*pi;...
               242   347   0/180*pi;...
               211   244  100/180*pi];
    if (pigMode == 2)
        meanPos = repmat([242, 347, 0/180*pi], 3, 1);
    end
  end

  % Model
  load(sprintf('%shog%dmark%d', resDir, mod(pigMode-1,2)+1, marking));
  likeParams.hogParams = hogParams;
  likeParams.cropRotatedIm = cropRotatedIm;

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % parameters for state space definition

  dim = 3;
  numsamples = 5000;


  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % params for initialization

  initMean = [meanPos(1,:); meanPos(3,:)];
  initCovar= diag([5^2, 5^2,(10/180*pi)^2]);

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % parameters for dynamics

  % covariance matrix of dynamical diffusion
  % YOU NEED TO CHANGE THIS
  % dynamics.dynCovar =  diag([1^2 1^2 (1/180*pi)^2]);
  
  if (pigMode == 1 || pigMode == 3)
    if (marking == 0) % Grey/Back
        dynamics.dynCovar = diag([4^2 4^2 (50/180*pi)^2]);
    else
        if (marking == 1) % Grey/Head
          dynamics.dynCovar = diag([4^2 4^2 (60/180*pi)^2]);  
        end
    end
  else
      if (pigMode == 2)
        if (marking == 0) % Red/Back
            dynamics.dynCovar = diag([3.5^2 3.5^2 (30/180*pi)^2]);
        else
            if (marking == 1) % Red/Head
                dynamics.dynCovar = diag([4^2 4^2 (40/180*pi)^2]);  
            end
        end
      end
  end

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Allocate and initialize sample sets to approximate distributions
  numsamples1 = round(numsamples/2);
  
  sampDist.samples = randn(numsamples, dim)*chol(initCovar) +  ...
      [repmat(initMean(1,:),numsamples1,1);
       repmat(initMean(2,:), numsamples - numsamples1, 1)];
  sampDist.weights = ones(numsamples,1) / numsamples;


  %% Loop over the image filenames in fn and go nuts
  frameNum0 = seqParams.frameStart;
  for frameNum = frameNum0:seqParams.frameSkip:seqParams.frameEnd
    % Read and display the frameNum^th image
    im = pgmRead([seqParams.frameRoot seqParams.framePrefix ...
                  num2strPad(frameNum, 4) '.pgm']);

    figure(1); clf; showIm(im);
    title(sprintf('Frame %d', frameNum));
    pause(0.1);
    
    if safetyNet  
      % Use coarse track results to get rough localization.
      % Crop rough neighbourhood of head region.
      bodyMn = coarseRes{frameNum}.mn(pigMode,:);
      U = coarseRes{frameNum}.U(:,:,pigMode);
      xc0 = round(bodyMn + U(:,1)' * mn2HeadDist);
    else
      % Use previous mean pose to get rough localization.
      meanPose = sampleMean(sampDist);
      xc0 = round(meanPose(1:2));
    end
    
    xRange = xc0(1) + ((-sqrCropRad):sqrCropRad);
    yRange = xc0(2) + ((-sqrCropRad):sqrCropRad);
    xRange = xRange(xRange >=1 & xRange <= size(im,2));
    yRange = yRange(yRange >=1 & yRange <= size(im,1));

    imCrop = im(yRange, xRange);
    figure(1); hold on;
    plot([xRange(1) xRange(1) xRange(end) xRange(end) xRange(1)],...
         [yRange(1) yRange(end) yRange(end) yRange(1) yRange(1)], 'g');
    
    % Display image and results
    figure(2); clf; showIm(imCrop);
    title('Cropped search region');
    
    % draw samples from posterior at time t-1
    sampDist1 = resampleRandom(sampDist);
    % sampleMean(sampDist1)
    
    % propagate samples through dynamics and resample
    sampDist2 = propagateSamples(sampDist1, dynamics);
    thetas = sampDist2.samples(:,3);
    idx = thetas > pi;
    while any(idx)
      thetas = thetas - 2*pi;
      idx = thetas > pi;
    end
    idx = thetas <= -pi;
    while any(idx)
      thetas = thetas + 2*pi;
      idx = thetas <= -pi;
    end
    
    % reweight the samples using the likelihood based importanc weights
    sampDist2.samples = sampDist2.samples - ...
        repmat([xRange(1)-1, yRange(1)-1, 0], numsamples, 1);
    sampDist3 = updateWeights(sampDist2, imCrop, likeParams);
    sampDist3.samples = sampDist3.samples + ...
        repmat([xRange(1)-1, yRange(1)-1, 0], numsamples, 1);
    
    sampleMean(sampDist3)
    

    % display the samples weighted by the likelihood
    figure(3); clf;
    showIm(log(displayDist(sampDist3,size(im))+1/numsamples));
    title('Log Filtering Distribution (position only)');

    meanPosterior = sampleMean(sampDist3) - [xRange(1), yRange(1), 0];
    figure(2); hold on;
    h = plot(meanPosterior(1), meanPosterior(2), 'or', 'MarkerSize', 4, ...
             'MarkerFaceColor', 'r');
    thetas = sampDist3.samples(:,3);
    r = mean([cos(thetas) sin(thetas)], 1);
    theta = atan2(r(2), r(1));
    tang = [cos(theta), sin(theta)];
    tang = 10*tang;
    h = plot(meanPosterior(1)+[0,tang(1)], meanPosterior(2)+[0,tang(2)],...
             'r', 'LineWidth', 2);
    
    % Display orientation distribution
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
    figure(4); clf;
    plot(thetaBin, thetaHist/sum(thetaHist), 'LineWidth',2);
    ax = axis;
    axis([0 2*pi 0, ax(4)]);
    set(gca, 'FontSize', 14);
    xlabel('Particle angle: \theta');
    ylabel('Relative frequency');
    title('Orientation distribution');

    sampDist = sampDist3;

    save(sprintf('results/%s%dmark%dfrm%s.mat', outPrefix, pigMode, ...
          marking, num2strPad(frameNum, seqParams.padNum)), 'sampDist');
    pause(0.1);  %% ensure plots are drawn
  end

end