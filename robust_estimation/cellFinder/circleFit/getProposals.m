% Alex Hong, 2014
function [circles] = getProposals(normals, p, numGuesses)
  % [circles] = getProposals(normals, p, numGuesses)
  % Attempt to produce up to numGuesses circle proposals from
  % the edgel data p and normals.  For typical data sets
  % we will be able to produce numGuesses proposals.  However,
  % on some datasets, say with only a few edgels, we may not be
  % able to generate any proposals.  In this case size(circles,1)
  % can be zero.
  % Input:
  %  normals - N x 2 edgel normals
  %  p         N x 2 edgel positions
  %  numGuesses - attempt to propose this number of circles.
  % Return:
  %   circles a P x 3 array, each row contains [x0(1) x0(2) r]
  %           with 0 <= P <= numGuesses.
  
  circles = zeros(0, 3); % proposals
  N = size(normals, 1);
  r_ratioThreshold = 1.3;
  
  % Bounding box
  xMin = min(p(:, 1));
  xMax = max(p(:, 1));
  yMin = min(p(:, 2));
  yMax = max(p(:, 2));
  
  % Maximum possible radius
  %r_maxThres = sqrt((xMax - xMin)^2 + (yMax - yMin)^2) / 2;
  r_maxThres = 35;
  
  % Check
  if (N < 2)
    fprintf('Error: Not enough edgels');
    return;
  end
  
  % Pair each edgel with another edgel and find their intersection
  for n = 1:numGuesses
    while(true)
         % get edgel indices randomly, e_i(1) != e_i(2)
         e_i = randsample(1:N, 2);
         p1 = e_i(1);
         p2 = e_i(2);
      
          % calculate slopes (rise / run)
          m(1) = normals(p1, 2) / normals(p2, 1);
          m(2) = normals(p2, 2) / normals(p2, 1);
          % check if normals are parallel
          if(m(1) ~= m(2))
             b(1) = p(p1, 2) - m(1) * p(p1, 1);
             b(2) = p(p2, 2) - m(2) * p(p2, 1);
             % find intersection of edgels
             x_intersect = (b(2) - b(1)) / (m(1) - m(2));
             y_intersect = m(1) * x_intersect + b(1);
             % check if the intersection point lies in the direction of
             % the two normals
             if ((normals(p1, 1) >= 0) == (x_intersect >= p(p1, 1)) && ...
                 (normals(p1, 2) >= 0) == (y_intersect >= p(p1, 2)) && ...
                 (normals(p2, 1) >= 0) == (x_intersect >= p(p2, 1)) && ...
                 (normals(p2, 2) >= 0) == (y_intersect >= p(p2, 2)))
                
                % find the two radii
                r(1) = pdist([p(p1, 1), p(p1, 2); x_intersect, y_intersect]);
                r(2) = pdist([p(p2, 1), p(p2, 2); x_intersect, y_intersect]);
                
                % checks if the calculated radii are near each other
                if (r(1) / r(2) <= r_ratioThreshold && ...
                    r(2) / r(1) <= r_ratioThreshold && ...
                    r(1) < r_maxThres && r(2) < r_maxThres )
                    
                    r = mean([r(1) r(2)]);
                    if (r < 25)
                        circles(n, :) = [x_intersect, y_intersect, r];
                        break;
                    end
                end
             end
          end
      end 
  end
end