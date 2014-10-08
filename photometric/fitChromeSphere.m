% Alex Hong (c) 2014
function [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % [L] = fitChromeSphere(chromeDir, nDir, chatty)
  % Input:
  %  chromeDir (string) -- directory containing chrome images.
  %  nDir -- number of different light source images.
  %  chatty -- true to show results images. 
  % Return:
  %  L is a 3 x nDir image of light source directions.

  % Since we are looking down the z-axis, the direction
  % of the light source from the surface should have
  % a negative z-component, i.e., the light sources
  % are behind the camera.
    
  if ~exist('chatty', 'var')
    chatty = false;
  end
    
  mask = ppmRead([chromeDir, 'chrome.mask.ppm']);
  mask = mask(:,:,1) / 255.0;

  for n=1:nDir
    fname = [chromeDir,'chrome.',num2str(n-1),'.ppm'];
    im = ppmRead(fname);
    imData(:,:,n) = im(:,:,1);           % red channel
  end

  % Light source direction
  L = zeros(3, nDir);
  
  % Center of sphere
  % calculate the center of mass of sphere with respect to
  % the origin (x, y) = (0, 0)
  P_center = [sum((1:size(mask, 2)) * mask') / sum(sum(mask));
              sum((1:size(mask, 1)) * mask) / sum(sum(mask));];
             
  % Radius of sphere
  % average of maximum height and maximum width of the sphere
  radius = (max(sum(mask')) + max(sum(mask))) / 4;
  
  % For each light source image
  for i = 1:nDir
     % Brightest spot within sphere
     img = (imData(:,:,i) .* mask).^2;
     P_brightest = [sum((1:size(mask, 2)) * img') / sum(sum(img));
                    sum((1:size(mask, 1)) * img) / sum(sum(img));];
                
     % Surface normal of brightest spot on sphere
     n_x = P_brightest(1) - P_center(1);
     n_y = P_brightest(2) - P_center(2);
     normal = [n_x;
               n_y;
               -sqrt(radius^2 - n_x^2 - n_y^2);];
     normal = normal / norm(normal); % normalized surface normal
     
     % Light source direction
     d_e = [0;
            0;
           -1;]; % camera direction
     
     L(:,i) = 2 * dot(d_e, normal) * normal - d_e; % Light source direction
     % Show results
     if (chatty)
         figure(1); clf;
         showIm(imData(:,:,n) .* mask);
     end 
  end
  return;

