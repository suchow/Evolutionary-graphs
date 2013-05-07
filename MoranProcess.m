function history = MoranProcess(N,K,numSteps,graphType)
  
  % default graph is complete
  if(nargin < 4)
    graphType = 'Complete';
  end
  
  % initialize the population
  population = MakePopulation(N,graphType{1},graphType{2:end});
  population = Endow(population,K,true);
  population = Endow(population,K,false,'random');
  
  % run the process
  history = cell(1,numSteps);
  for stepIndex = 1:(numSteps+1)
    
    % select an individual with probability proportional to fitness
    i = Select(population,0);
    
    % look at the selected individual's outgoing edges, and choose
    % another individual who this individual's offspring will replace
    j = randp(full(population.graph(i,:)));

    % record the history
    history{stepIndex} = population; 
    history{stepIndex}.death = j;
    history{stepIndex}.reproduce = i;
    
    % replace the individual
    population.tags(j) = population.tags(i);
    population.memories(j,:) = population.memories(i,:);
  end
end
