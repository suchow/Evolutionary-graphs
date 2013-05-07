function population = Endow(population,K,isRandomAllocation,virality)

    if(nargin < 3)
      isRandomAllocation = true;
    end
    
    if(nargin < 4)
      virality = 'random';
    end

    N = size(population.graph,1);

    switch virality

      case 'random'
        % pick types (i.e., which of the K events each individual observed)
        if(isRandomAllocation)
          population.tags = randi(K,N,1);
        else
          population.tags = sort(1+mod(0:N-1, K))';
        end
    
      case 'infectious'
        % place the first K observations randomly
        population = placeK(population);
      
        % infect
        while(any(population.tags == 0))
          j = randsample(find(population.tags));
          q = find(mnrnd(1,full(population.graph(j,:))));
          population.tags(q) = population.tags(j);
        end
        population.tags = population.tags';
      
      case 'nearest'
        % place the first K observations randomly
        population = placeK(population);
        seedIndices = find(population.tags);
        
        % compute path lengths
        pathlengths = pathlength(population.graph);
        
        % for each of the remaining N-K observations, assign by copying their 
        % nearest neighbor on the graph.
        toFill = find(population.tags == 0);
        for i = 1:length(toFill)
          distanceToSeeds = pathlengths(toFill(i),seedIndices); % compute distance to seeds
          [~,nearestSeedIndex] = mintb(distanceToSeeds);
          population.tags(toFill(i)) = population.tags(seedIndices(nearestSeedIndex));
        end
        population.tags = population.tags';
    end
    
    % pick memories. memories are made up of two numbers: an event assignment
    % and a unique identifier.
    population.memories = [population.tags, (1:N)'];
    
    
    function pop = placeK(pop)
      % place the first K observations randomly
      ord = shuffle(1:N);
      pop.tags = zeros(1,N);
      for k = 1:K
        pop.tags(ord(k)) = k;
      end
    end
end