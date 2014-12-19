% Alex Hong, 2014
% Script file:  Extract a hog model from one image.
%  Mouse in the pigeon head center and orientation (or use
%  the given points).
%  Rotate image.  Extract HOG model centered on moused in point.
clear;
close all;
fclose all;
clc;

run '../../matlab/startup.m'

onCDF = false;
if onCDF
  seqRoot = '/u/kyros/public/csc2503/pigeon/';
else
  % Change this to your sequence install directory
  % if you are not on CDF.
  seqRoot = '../pigeon/';
end

marking = 0  % 0 denotes back marking, 1 denotes head
useMouse = false;

% pigMode 1/2/3 denotes :
%    1  the grey pigeon, top right in frame 30,
%    2  the red pigeon
%    3  the cloned grey pigeon, middle left in frame 30.
  
for pigMode = 1:2  % Use 1 or 2 or 1:2
  pigMode
  frameNum = 30;
  im = pgmRead([seqRoot 'pigClone' num2strPad(frameNum, 4) '.pgm']);
  figure(1); clf; showIm(im);
  if (pigMode == 1)
    % Pigeon 1: (Blue)
    figure(1);
    % Blow up image before mousing in points.
    if useMouse
      pause;
      xc = round(ginput(1)); 
      beak = ginput(1); 
    else
      if marking == 1  % head 
        xc = [470   115];
        beak = [449   109];
      elseif marking == 0  % back
        xc = [531   146];
        beak = [480   112];
      end
    end 
  else
    figure(1);
    % Blow up image before mousing in points.
    if useMouse
      pause;
      xc = round(ginput(1)); 
      beak = ginput(1); 
    else
      if marking == 1  % head 
        xc = [243  348];
        beak = [261   350];
      elseif marking == 0  % back
        xc = [158   322];
        beak = [210   341];
      end

    end 
    
  end

  % Size of warped cropped image is (2*rad+1) square.
  rad = 100;
  xcw = (rad +1) * ones(1,2);


  % Build similarity warp from the canonical image to this image.
  A = eye(3);
  ax = (beak-xc);
  ax = ax/norm(ax);

  % R e1 = ax
  R = [ax' [-ax(2); ax(1)]];
  A(1:2, 1:2) = R;

  % A * [xcw, 1]' = [xc, 1]'
  A(1:2,3) = xc' - R * xcw';

  % Do the warp
  [imB, B] = imWarpAffine(im, A, 1);
  imB(isnan(imB)) = 0;

  % Crop the result.
  imCrop = imB(1:(2*rad+1), 1:(2*rad+1));

  figure(2); clf; showIm(imCrop);

  z = (rad+1) * ones(2,1);
  hogParams.maxAmp = 25;
  hogParams.minAmp = 6;
  hogParams.nAmp = 3;
  hogParams.nTheta = 10;
  hogParams.nxCell = 3;
  hogParams.nyCell = 3;
  hogParams.sx = 15;
  hogParams.sy = 10;
  hogParams.lambda = 10;
  wCell = ones(hogParams.nyCell, hogParams.nxCell);
  hogParams.wCell = wCell;

  histModel = getAlignedHOG(z, imCrop, hogParams, true);

  histModel = reshape(histModel, [hogParams.nyCell, hogParams.nxCell, ...
                    hogParams.nAmp * hogParams.nTheta]);

  hogParams.histModel = histModel;
  hogParams.pigMode = pigMode;

  if true % Save the HOG model.
    if marking == 0
      hogDir = 'back/';
    elseif marking == 1
      hogDir = 'head/';
    end
    save(sprintf('%shog%dmark%d', hogDir, pigMode, marking), '-V6', 'hogParams');
  end
end
