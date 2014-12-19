function [hist] = getAlignedHOG(z, imAligned, params, display)   
  if nargin < 4
    display = false;
  end
  
  CHATTY = true;
  if display | CHATTY
    figure(201); clf; displayImage(imAligned, [0 255]); 
    hold on;
    plot(z(1,:), z(2,:), '*r');
    title('Aligned HOG centers');
    pause(0.1);
  end
  
  nPose = size(z,2);
  szIm = size(imAligned);
  
  % Bin boundaries in amplitude
  maxAmp = params.maxAmp;
  minAmp = params.minAmp;
  nAmp = params.nAmp;    
  minLAmp = log(minAmp);
  maxLAmp = log(maxAmp);
  sLAmp = (maxLAmp - minLAmp)/nAmp;
  binLAmp = (0:nAmp) * sLAmp + minLAmp;
  
  % Shape and number of HOG cells
  sx = params.sx;         
  sy = params.sy;        
  nxCell = params.nxCell;
  nyCell = params.nyCell;
  

  % Direction image in range (-1,1]
  nTheta = params.nTheta;  % 10
  sDir = 2/nTheta;
  binDir = ((1:nTheta)-1) * sDir;
  binDir(binDir>1) = binDir(binDir>1) - 2;

  % Initialize histogram
  hist = zeros([nxCell, nyCell, nTheta, nAmp, nPose]);

  % Gradient image
  sigmaGrad = 1.0;
  sigmaSqr = sigmaGrad*sigmaGrad;
  gGradSize = 2 * round(3.0 * sigmaGrad) + 1;
  x = [1:gGradSize] - round((gGradSize+1)/2);
  gGrad = exp(- x .* x / (2.0*sigmaSqr));
  gGrad = gGrad/ sum(gGrad(:));
  gxGrad = - (x / sigmaSqr) .* gGrad;

  gradIm = zeros(size(imAligned,1), size(imAligned,2), 2);

  gradIm(:,:,1) = rconv2sep(imAligned, gxGrad, gGrad);
  gradIm(:,:,2) = rconv2sep(imAligned, gGrad, gxGrad);

  % Direction of gradient image
  dirGrad = atan2(gradIm(:,:,2), gradIm(:,:,1))/pi;
  
  % Log of norm of gradient image
  logNrmGradIm = log(sum(gradIm.^2, 3)+eps)/2.0;
  
  % Saturate the log amplitude at the maximum binLAmp(end)
  logNrmGradIm = min(binLAmp(end), logNrmGradIm);
  
  % Threshold small gradients to zero.
  idx = logNrmGradIm <= binLAmp(1);
  if (any(any(idx)))
    logNrmGradIm(idx) = 0;
  end
  
  % Build the integral images for each amplitude range
  % and gradient orientation range.
  for ka = 2:length(binLAmp)
    
    % Find all points with gradient log amplitude in the
    % range (binLAmp(ka-1), binLAmp(ka)] 
    lgBin = zeros(size(logNrmGradIm));
    idx = (logNrmGradIm >binLAmp(ka-1)) & (logNrmGradIm <= binLAmp(ka));
    if any(any(idx))
      lgBin(idx) = logNrmGradIm(idx);
    end
    
    if false & CHATTY
      figure(203); clf; showIm(lgBin); pause(0.1);
    end
    
    % Find gradients within the range of a HOG bin.
    for kd = 1:nTheta
      
      % Use the directions in [binDir(kd)-sDir/2, binDir(kd)+sDir/2)
      % (Recall all gradient directions have been divided by pi.)
      thetaGrad = zeros(size(logNrmGradIm));
      diffDir = dirGrad - binDir(kd);
      idx = diffDir >= 1;
      while (any(any(idx)))
        diffDir(idx) = diffDir(idx) - 2;
        idx = diffDir >= 1;
      end
      idx = diffDir < -1;
      while (any(any(idx)))
        diffDir(idx) = diffDir(idx) + 2;
        idx = diffDir < -1;
      end
      idx = (lgBin > 0) & (diffDir >= -sDir/2.0) & (diffDir < sDir/2.0);
      if (any(any(idx)))
        thetaGrad(idx) = lgBin(idx);
      end
      
      if false & CHATTY
        figure(204); clf; showIm(thetaGrad);
        title(sprintf('HOG bin: a (%3.1f, %3.1f], tb [%3.2f, %3.2f]',...
           exp(binLAmp((ka-1):ka)), pi*(binDir(kd) + sDir/2.0 *[-1, 1])))  
        pause(0.1);
      end
      
      % Pad the gradient data with a row and column of zeros.
      padThetaBin = zeros(size(thetaGrad)+1);
      padThetaBin(2:end, 2:end) = thetaGrad;
      
      % Get the integral image
      intIm = integralIm(padThetaBin);
      szIntIm = size(intIm);
      intIm = intIm(:);
      
      % Extract the histogram bins for this grad amp and orientation.
      xcg = round((1 + [nxCell * sx, nyCell * sy])/2.0);
      for iy = 1:nyCell
        for ix = 1:nxCell
          
          % Histogram spatial extent [x0, x1) x [y0, y1)
          x0 = round(z(1,:) + (1+(ix-1)*sx - xcg(1)));
          y0 = round(z(2,:) + (1+(iy-1)*sy - xcg(2)));
          x1 = x0 + sx;
          y1 = y0 + sy;
          
          % Clip to boundaries of image
          x0 = max(x0, 1);
          y0 = max(y0, 1);
          x1 = max(x0, x1);
          y1 = max(y0, y1);
          
          x1 = min(x1, szIm(2));
          y1 = min(y1, szIm(1));
          x0 = min(x0, x1);
          y0 = min(y0, y1);
          
          if  display & ka == 2 & kd == 1
            for k=1:length(x0)
              figure(201); hold on;
              plot([x0(k) x1(k)], [y0(k) y0(k)],'r');
              plot([x0(k) x1(k)], [y1(k) y1(k)],'r');
              plot([x0(k) x0(k)], [y0(k) y1(k)],'r');
              plot([x1(k) x1(k)], [y0(k) y1(k)],'r');
            end
            % print(figure(201), '-depsc', 'pigModel');
          end
 
          % Build value within spatial extent.
          val = intIm( y1 + szIntIm(1)* (x1-1));
          val = val - intIm(y0 + szIntIm(1) * (x1 - 1));
          val = val - intIm(y1 + szIntIm(1) * (x0-1));
          val = val + intIm(y0 + szIntIm(1) * (x0 - 1));
          hist(iy, ix, kd, ka-1, :) = val;
        end % over ix
      end  % over iy
    end % over kd
  end % over ka

  hist = max(hist, 0.0);
  return;
  
