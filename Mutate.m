function population = Mutate(population,p)
  % shall we?
  if(rand < p)
    i = randi(size(population.graph,1));
    population.anyMutations = true;
    population.memories(i,2) = population.maxMutation + 1;
    population.maxMutation = population.maxMutation + 1;
  end
end