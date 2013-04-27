% tags of -1 are inaccessible, 1 through N are normal, > N+1 are noise
function s = MemoryStatistics(population)
  
  % compute number remembered
  s.numRemembered = length(unique(population.tags(population.tags ~= -1)));
  isNoise = (population.memories(:,2) > population.N | ...
             population.memories(:,2) == -1);
  s.numSurvivingLineages = length(unique(population.memories(~isNoise,2)));
end