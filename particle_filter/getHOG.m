function [hist] = getHOG(pose, im, params) 
% [hist, hogParams] = getHOG(pose, im, params) 
% Build histogram of gradients for a set of poses.
% Input: 
%   pose 3 x np matrix with columns (xc, yc, theta)' containing
%        the proposed head coordinates in the coordinates of the provided
%        image im.
%   im the ny x nx grey-level image to search.
%   params.hogParams  the HOG model params and target histogram
%   params.cropRotatedIm true [Default] to crop rotated iamge.
% Returns:
%  hist (nAmp*nTheta) x np HOG histograms
%  
  hogPar = params.hogParams;
  if params.cropRotatedIm
    cropRotIm = 1;
  else
    cropRotIm = 0;
  end
  
  nThetaCrop = 24;
  dTheta = 2*pi/24.0;
  thetaCrop= dTheta * (0:(nThetaCrop-1));
  
  xc0 = ([size(im,2), size(im,1)]+1)/2.0;
  sqrCropRad = floor(min(size(im))/2.0);
  
  nPose = size(pose,2);
  
  % Find the closest discrete orientation for each pose
  theta = pose(3,:);
  thetaCrop = thetaCrop(:);
  nTheta = length(thetaCrop);
  dTheta = repmat(theta, nTheta, 1) - repmat(thetaCrop, 1, nPose);
  
  % Do mod 2*pi on the difference
  dTheta = dTheta(:);
  idx  = (dTheta > pi);
  while (any(idx))
    dTheta(idx) = dTheta(idx) - 2*pi;
    idx = (dTheta > pi);
  end
  idx  = (dTheta <= -pi);
  while (any(idx))
    dTheta(idx) = dTheta(idx) + 2*pi;
    idx  = (dTheta <= -pi);
  end
  dTheta = reshape(dTheta, nTheta, nPose);
  [tmp, kTheta] = min(abs(dTheta), [], 1);

  % Initialize results to zero
  szHist = prod([hogPar.nxCell hogPar.nyCell hogPar.nTheta hogPar.nAmp]);
  hist = zeros(szHist, nPose);
   
  % Prepare for looping over rotations
  szIm = size(im); 
  
  % Loop over discrete image rotations
  for kt = 1:nThetaCrop
    % Find poses rounded to this image rotation
    idx = kTheta == kt;
    if (any(idx))
      % Do the image rotation
      thetaC = thetaCrop(kt);
      % Rotate image by the angle -thetaC (counter-clockwise).
      % Build the image transform.
      ct = cos(thetaC);
      st = sin(-thetaC);
      R = [ct st; -st ct];
      A = eye(3); A(1:2,1:2) = R; 
      % Set the center of the final warped cropped images
      % Make A*[sqrCropRad, sqrCropRad, 1]' = [xc0', 1]'
      A(1:2,3) = xc0(:) - R*[sqrCropRad ; sqrCropRad];
      
      % Do the warp
      [imB, B] = imWarpAffine(im, A, cropRotIm);
      imB(isnan(imB)) = 0;
      
      % Crop the result.
      imCrop = imB;
      
      % Map the pose points according to the rotation.
      z = inv(B) * [pose(1:2, idx); ones(1, sum(idx))];
      
      % Get the HOG's for all poses closest to this image rotation.
      qw =  getAlignedHOG(z, imCrop, hogPar);
      hist(:,idx) = reshape(qw, szHist, sum(idx));
    end
  end
    
  return;