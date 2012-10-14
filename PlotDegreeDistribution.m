function PlotDegreeDistribution(graphType, N, varargin)
  numIterations = 1;
  if(nargin < 2)
    N = 1000;
  end
  data = zeros(numIterations,N);
  parfor i = 1:numIterations
    fprintf('%d\n', i)
    w = MakeAdjacencyMatrix(graphType,N,varargin{:})
    data(i,:) = histc(full(sum(w>0)), [1:N]);
  end
  loglog([1:N], mean(data,1)./N, '.k')
  makepalettable
  xlabel('Degree')
  ylabel('Probability')