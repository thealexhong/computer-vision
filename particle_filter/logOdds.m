function [lgOdds] = logOdds(resp, lambda, mu, scl)
  
  qw = exp(-lambda * (resp - mu));
  lgOdds = scl * (1- qw)./(1 + qw);
  return;