% frequency dependent selection with parameter tau that determines the
% strength of selection. tau is not yet implemented.
function i = Select(population, policy)

  if(nargin < 2)
    i = randsample(find(~whererare(population.tags,1)));
  else
    i = policy(population);
  end
end
