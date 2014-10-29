% Alex Hong, 2014
function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr, normals, sigmaGM)
%
% function [x0, r, w, maxW] = fitCircleRobust(pts, initx0, initr,
%                                  normals, sigmaGM)
%
%  minimize sum_i  rho( a'x_i + b + x_i'x_i, sigma)
%  w.r.t. a and b, where a = -2x_0 and b = x_0'x_0 - r^2
% Input:
%  pts: an Nx2 matrix containing 2D x_i, i=1...N
%  initx0: 2-vector initial guess for centre
%  initr: initial guess for radius 
%  normals: an Nx2 matrix of corresponding 2D edge normal directions nrml_i  
% Output
%  x0 - 2x1 fitted circle's center position
%  r  - fitted circle's radius
%  w  - N x 1 robust weights for edgels.
%  maxW  - maximum possible robust weight (i.e. for zero error)
%          This may be used by the calling code to scale the weights 
%          to be between 0 and 1.
  
  x0 = initx0(:)';  % Make it a row vector.

  K = size(pts, 1);
  x_c = initx0; % initial circle centre guess
  r_k = initr; % initial radius guess
  x_k = pts; % edgels K x 2
  param_sumSQR0 = 0; % convergence variable
  
  % Loop until convergence
  %while(true)
  for iter = 1:50
    % Calculate errors of edgels
    x_cSQR = sum(x_c .^ 2); % 1 x 1
    x_cSQR = repmat(x_cSQR, 1, K); % 1 x K
    
    x_kSQR = sum((x_k .^ 2)')'; % K x 1
    
    
    r_k = repmat(r_k, 1, K); % 1 x K
    
    a = -2 * x_c; % 2 x 1
    b = x_cSQR - r_k .^ 2; % 1 x K
    b =  b';
    e_k = x_k * a + b + x_kSQR; % K x 1
    
    % Calculate weight of edgels
    e_kSQR = e_k .^ 2; % K x 1
    w = (sigmaGM./(sigmaGM^2 + e_kSQR)).^2; % K x 1
    W = diag(w); % K x K
    maxW = max(w);
    
    % Construct X matrix, defined in report
    X = ones(K, 3);
    X(:, [1, 2]) = pts; % K x 3
    
    % Solve for param, equation derived in report
    param = (X' * W * X) \ (-X' * W * x_kSQR);
    a = param([1 2]);
    b = param(3);
    x_c = a ./ (-2);
    x_cSQR = sum(x_c .^ 2);
    r_k = sqrt(x_cSQR - b);
    
    % check for convergence
    param_sumSQR = sum(param);
    if abs(param_sumSQR - param_sumSQR0) < 0.05
        break;
    end
    param_sumSQR0 = param_sumSQR;
  end
  % outputs
  x0 = x_c;
  r = r_k;
end

