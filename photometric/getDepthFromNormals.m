% Alex Hong (c) 2014
function [depth] = getDepthFromNormals(n, mask)
  % [depth] = getDepthFromNormals(n, mask)
  %
  % Input:
  %    n is an [N, M, 3] matrix of surface normals (or zeros
  %      for no available normal).
  %    mask logical [N,M] matrix which is true for pixels
  %      at which the object is present.
  % Output
  %    depth an [N,M] matrix providing depths which are
  %          orthogonal to the normals n (in the least
  %          squares sense).
  %
  imsize = size(mask);
  nPix = size(mask(:), 1);
  
  % Construct the A variable
  nz = reshape(n(:,:,3), [nPix, 1]);
  Ay = spdiags([-nz nz],[0 1],nPix,nPix);
  Ax = spdiags([-nz nz],[0 imsize(1)],nPix,nPix);
  A = [Ax;
       Ay;];
  
  
  % Construct the v variable
  nxy = [reshape(-n(:,:,1), [nPix, 1]);
         reshape(-n(:,:,2), [nPix, 1]);];
  v = nxy;
  
  
  depth = A\v;
  depth = reshape(depth, imsize);

  
end


