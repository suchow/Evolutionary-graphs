function population = Endow(population,K,isRandomAllocation,virality)

    if(nargin < 3)
      isRandomAllocation = true;
    end
    
    if(nargin < 4)
      virality = 'none';
    end

    N = size(population.graph,1);

    switch virality

      case 'none'
        % pick types (i.e., which of the K events each individual observed)
        if(isRandomAllocation)
          population.tags = randi(K,N,1);
        else
          population.tags = sort(1+mod(0:N-1, K))';
        end
    
      case 'spread'
        % place the first K observations
        ord = shuffle(1:N);
        population.tags = zeros(1,N);
        for i = 1:K
          population.tags(ord(i)) = i;
        end
      
        % infect
        while(any(population.tags == 0))
          j = randsample(find(population.tags));
          q = find(mnrnd(1,full(population.graph(j,:))));
          population.tags(q) = population.tags(j);
        end
        population.tags = population.tags';
    end
    
    % pick memories. memories are made up of two numbers: an event assignment
    % and a unique identifier.
    population.memories = [population.tags, [1:N]'];
end