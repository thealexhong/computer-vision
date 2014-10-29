% Alex Hong, 2014
function [circle] = bestProposal(circles, sigmaGM, normals, p)
% [circle] = bestProposal(circles, sigmaGM, normals, p)
% Chose the best circle from the proposed circles.
% Input
%  circles K x 3 matrix each row containing [x0(1), x0(2), r] for 
%          the center and radius of a proposed circle.
%  sigmaGM - the robust scale parameter 
%  normals - P x 2 edgel normal data
%  p       - P x 2 edgel position data.
% Output
%  circle  1 x 3  best circle params [x0(1), x0(2), r] from the proposals
  
  K = size(circles, 1);
  P = size(normals, 1);
  
  % Error is defined as e_k = ||x_k - x_c||^2 - r^2
  
  % get x_c and x_c^T*x_c for P edgel positions
  x_c = circles(:,[1 2]); % K x 2
  x_cSQR = sum(x_c.^2,2); % K x 1
  x_cSQR = repmat(x_cSQR, 1, P); % K x P
  
  % get x_k and x_k^T*x_k for K proposals
  x_k = p; % P x 2
  x_kSQR = sum(x_k.^2, 2); % P x 1
  x_kSQR = repmat(x_kSQR, 1, K); % P x K

  % get R for P edgel positions
  r = circles(:, 3); % K x 1
  r = repmat(r, 1, P); % K x P
  
  % Parameterization of a = -2 * x_c; b = x_c^T * x_c - r^2
  a = -2 * x_c; % K x 2
  b = x_cSQR - r.^2; % K x P
  
  % Error
  e_k = x_k * a' + b' + x_kSQR; % P x K
  
  % pick the circle with the least error
  [~, i] = min(sum(e_k .^ 2));
  circle = circles(i,:);
end