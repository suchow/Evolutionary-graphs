% min with a tie breaker among equal options
function [y,i] = mintb(x)
  [y,i] = min(x);
  matchIndices = find(x==y);
  i = randsample(matchIndices);
end