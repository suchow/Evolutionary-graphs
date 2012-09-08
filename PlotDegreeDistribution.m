function PlotDegreeDistribution(graphType, N, varargin)
  numIterations = 10;
  if(nargin < 2)
    N = 1000;
  end
  data = zeros(numIterations,N);
  parfor i = 1:numIterations
      w = MakeAdjacencyMatrix(graphType,N,varargin{:})
      data(i,:) = sort(full(sum(w>0)),'descend');
  end
  loglog([1:N], mean(data), '.k')
