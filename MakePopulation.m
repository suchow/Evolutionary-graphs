function population = MakePopulation(N,graphType,varargin)
    
    % default to a complete graph unless specified
    if(nargin < 2)
        graphType = 'Complete';
    end

    population.graphType = graphType;
    population.graph = MakeAdjacencyMatrix(graphType,N,varargin{:});
    population.N = size(population.graph,1);
    population.anyMutations = false;
    population.maxMutation = N;
end