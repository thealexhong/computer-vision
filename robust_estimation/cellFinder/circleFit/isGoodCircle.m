% Alex Hong, 2014
function [goodCircle] = isGoodCircle(x0, r, w, ...
                                     circleEstimates, nFound)
  % [goodCircle] = isGoodCircle(x0, r, w, normals, ...
  %                                  circleEstimates, nFound)
  % Decide if the circle with parameters x0 and r, with fitted
  % weights w, is to be added to the current set of circles.
  % Input:
  %  x0 2-vector for circle center
  %  r  radius of circle
  %  w  robust weights of edgels supporting this circle,
  %     weights have been scaled to have max possible value 1.
  %  circleEstimates C x 3 array of circle parameters, the first
  %     nFound rows of which specifies the [x0(1), x0(2) r] parameters
  %     for a previously fitted circle for this data set.
  %  nFound the number of previously fitted circles stored in circleEstimates
  % Output:
  %  goodCircle boolean, true if you wish to add this circle to the
  %             list of previously fitted circles.
  
  x0 = x0(:)';  % Decide, row or column
  
  wThreshold = 0.3;
  dThreshold = 0.25;
  
  % compare with threshold for 0 intersection
  if(nFound == 0)
      if (sum(w) > 2 * pi * r * wThreshold)
          goodCircle = true;
      else
          goodCircle = false;
      end
      return;
  end
  
  % Calculate circle intersections
  dist = zeros(nFound, 3);
  for i = 1:nFound
      d = pdist([x0 ; circleEstimates([1 2],i)']);
      dist(i, 1) = r + circleEstimates(3, i) - d;
      dist(i, 2) = circleEstimates(3, i); % radius of circleEstimates
  end
  
  crop = dist(:, 1) > 0; % intersections happen > 0
  dist = dist(crop, :);
  crop = dist(:, 1) < r * 0.6; % filter
  dist = dist(crop, :);
  dist(:, 3) = acos((r^2 + dist(:, 1).^2 - dist(:, 2).^2) ./ ...
                    (2 * dist(:, 1) .* dist(:, 2)));
  angle_diff = (2 * pi) - 2 * sum(dist(:, 3));
  % compare with threshold
  if (sum(w) > angle_diff * r * dThreshold)
      goodCircle = true;
  else
      goodCircle = false;
  end
end