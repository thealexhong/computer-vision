% Alex Hong (c) 2014
function [n, albedo] = fitReflectance(im, L)
  % [n, albedo] = fitReflectance(im, L)
  % 
  % Input:
  %   im - nPix x nDirChrome array of brightnesses,
  %   L  - 3 x nDirChrome array of light source directions.
  % Output:
  %   n - nPix x 3 array of surface normals, with n(k,1:3) = (nx, ny, nz)
  %       at the k-th pixel.
  %   albedo - nPix x 1 array of estimated albdedos
  
  g = (L * L')\(L * im');
  albedo = sqrt(sum(g.^2))';
  %n = g' ./ repmat(albedo, 1, 3);
  n = normr (g');
  
  return;


