function loglikelihood = pigLogLike(state, im, likeParams)
  
  hogParams = likeParams.hogParams;
  pigMode = hogParams.pigMode;
  nPose = size(state,2);
  histModel = hogParams.histModel;
  wCell = hogParams.wCell;
  
  hist = getHOG(state, im, likeParams);
  
  hist = reshape(hist, [hogParams.nyCell, hogParams.nxCell, ...
                    hogParams.nTheta*hogParams.nAmp, nPose]);

  sumX2 = zeros(nPose,1);
  cellResp = zeros(hogParams.nyCell, hogParams.nxCell, nPose);
  for ix = 1:hogParams.nxCell
    for iy = 1:hogParams.nyCell
      tmp = sum(hist(iy,ix,:,:).* ...
                repmat(histModel(iy, ix,:), [1,1,1,nPose]),3)./ ...
            sqrt((sum(hist(iy,ix,:,:).^2,3)+1)*...
                 (sum(histModel(iy,ix,:).^2)+1));
      cellResp(iy,ix,:) = tmp;
      sumX2 = sumX2 + reshape(wCell(iy, ix) * tmp, nPose, 1);
    end
  end

  resp = sumX2/sum(wCell(:));

  loglikelihood = logOdds(resp, hogParams.lambda, hogParams.mu, hogParams.scl);
 
  return;