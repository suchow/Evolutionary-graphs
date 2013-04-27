function history = ThresholdedMoranProcess(N,K,numSteps,graphType,threshold)
  
  % default graph is complete
  if(nargin < 4)
    graphType = 'Complete';
  end
  
  if(nargin < 5)
    threshold = 1;
  end
  
  % initialize the population
  population = MakePopulation(N,graphType{1},graphType{2:end});
  population = Endow(population,K,true,'spread');
  
  % run the process
  for stepIndex = 1:(numSteps+1)
    
    % select an individual with probability proportional to fitness
    i = Select(population,0,3);
    
    % look at the selected individual's outgoing edges, and choose
    % another individual who this individual's offspring will replace
    j = find(logical(mnrnd(1,full(population.graph(i,:)))));

    % record the history
    history{stepIndex} = population; 
    history{stepIndex}.death = j;
    history{stepIndex}.reproduce = i;
    
    % replace the individual
    population.tags(j) = population.tags(i);
    population.memories(j,:) = population.memories(i,:);
    
    % kill all memories represented by fewer than 'threshold' things
    r = whererare(population.tags, threshold);
    population.tags(r) = -1;
    population.memories(r,:) = repmat([-1,-1], sum(r), 1);
  end
end

% fast mnrnd 
function r = mnrnd(n,p)
  r = zeros(size(p));
  r(find(rand < cumsum(p),1)) = 1;
end