% Alex Hong, 2014

function [err, thres, parity, stumpResp, r] = trainStump(f, useAbs, X, y, ...
                                                  w, debug)
  %[err, thres, parity, stumpResp, r] = trainStump(f, useAbs, X, y, ...
  %                                                w, debug)
  % Input:
  %  f - N dimensional projection vector.
  %  useAbs - boolean, whether to use f(:)' * X or abs(f(:)' * X) as
  %     the feature.
  %  X  - N x K training matrix.
  %  y  - 1 x K label matrix, with y(k) = 1 indicating target, 0
  %       non-target.
  %  w  - 1 x K vector of weights.  Precondition: w(k) >= 0, sum_k w(k) > 0
  %  debug = 0, 1, or 2 for no info, histogram plot, more info,
  %         respectively.
  % Output
  %  err - (double) minimum weighted error of weak classifier (see (2) in
  %        the assignment handout).
  %  thres - (double) the optimum threshold theta_1 in (1) of the assignment.
  %  parity - (-1 or 1) the optimum parity theta_2 in (1) of the assignment.
  %  stumpResp - binary, 1 by K, optimum weak classifier response on X,
  %             i.e., the value of h(x,theta) in (1) of the assignment, 
  %             for the optimum parameters theta.
  %  r  - double 1 by K, the continuous valued response u(f(:)' * X)
  %       used by the optimum weak classifier (see (1) in the assignment).

  sumW = sum(w); % sum of all weights
  
  sumWNon = sum(w .* (y == 0)); % sum of weights for non-target
  
  r = f(:)' * X; % response
  if useAbs
    r = abs(r);
  end
  
  if debug & any(y > 0) & any(y == 0)
    % Plot histogram of (unweighted) responses.
    [hTarg, hxTarg] = histo(r(y>0), 101);
    hTarg = hTarg/(sum(hTarg) * (hxTarg(2) - hxTarg(1)));
    [hNon, hxNon] = histo(r(y==0), 101);
    hNon = hNon/(sum(hNon) * (hxNon(2) - hxNon(1)));
    figure(100); clf;
    plot(hxTarg, hTarg, 'g');
    hold on;
    plot(hxNon, hNon, 'r');
    xlabel('Response');
    ylabel('Probability');
    pause(0.1);
  end
  
  [sortR i] = sort(r); % sorts and grabs index of each response
  cumW1 = cumsum(w(i) .* y(i)); % cumulative sum of weighted targets
  cumW0 = cumsum(w(i) .* ~y(i)); % cumulative sum of weighted non-targets
  % for parity = 1:
  % error is sum of weight of non-targets under threshold + sum of weight
  % of targets over threshold
  errP = cumW0 - cumW1 + sumW - sumWNon;
  % for parity = -1:
  % error is sum of weight of targets under threshold + sum of weight
  % of non-targets over threshold
  errN = cumW1 - cumW0 + sumWNon;
  [errPmin errP_i] = min(errP);
  [errNmin errN_i] = min(errN);
  if (errPmin < errNmin) % find the minimum error
      parity = 1;
      err = errPmin / sumW;
      thres = r(i(errP_i));
      stumpResp = r <= thres;
  else
      parity = -1;
      err = errNmin / sumW;
      thres = r(i(errN_i));
      stumpResp = r > thres;
  end
  
  if debug
    % Plot threshold on histogram figure
    figure(100); hold on;
    ax = axis;
    plot([thres, thres], ax(3:4), 'k');
    title('Hist. of Weak Classifier, Targ/Non (g/r), Threshold (k)');
    pause(0.1);
  end

  return;

  
