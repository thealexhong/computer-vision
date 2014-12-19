%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File: tryEyeDetector.m
%  Matlab script file
%  Date: Dec, 2014
%

% Dependencies, Toolboxes:
%      iseToolbox/
% Author: Alex Hong, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

run '../../matlab/startup.m';
if ~exist('buildGaussFeat','file')
  addpath('./util');
end

which pgmRead

im = pgmRead('tryEye.pgm');

figure(1); clf;
showIm(im);
hold on;
szIm = size(im);
load('adaFit');

for nFaces = 1:6
    [x_i, y_i] = ginput(2);% select top left and bottom right
    x_i = floor(x_i);
    y_i = floor(y_i);
    width = ([-1 1] * x_i + 1);
    height = ([-1 1] * y_i + 1);
    szBox =  height * width; %N(height)xM(width)
    
    % patch coordinates relative to centers
    [x, y] = meshgrid((-9:10),(-12:12)); % 20x25 = 500 elements
    nElem = size(x,1) * size(x,2); % 500 elements
    % All patches' coordinates
    % x components of patch centres
    patchx = repmat((x_i(1):x_i(2)), height, 1); % NxM
    % y components of patch centres
    patchy = repmat(y_i(1):y_i(2), nElem, width); % 500 x N*M
    patchtotal = patchy + repmat(patchx(:)' * szIm(1), nElem, 1); % 500 x N*M
    % vector of patch coordinates
    patchcoord = x(:) * szIm(1) + y(:); % 500 x 1
    coords = patchtotal + repmat(patchcoord, 1, szBox); % 500 x N*M
    
    % Store patches' image intensities
    I = im(:);
    X = I(coords);
    % find eyes
    nPlot = size(featList, 2);
    resp = evalBoosted(featList, nPlot, X);
    eye = reshape(resp > 21.5, height, width);
    
    
    % draw bounding box and eye candidates on image
    [yEye, xEye] = find(eye);
    scatter(x_i(1) + xEye' - 1, y_i(1) + yEye' - 1)
    rectangle('Position',[x_i(1), y_i(1), width, height], ...
        'EdgeColor', 'g', 'LineWidth', 2)
end


