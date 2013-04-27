% takes in an adjacency matrix and returns graph stats
function stats = GraphStatistics(w,nargin)
  
  %
  % Local clustering coefficient
  %
  parfor i = 1:length(w) % find the neighbors
    neighborhood = find(w(i,:));
    numNeighbors = length(neighborhood);
    neighborhood(ismember(neighborhood, i)) = [];
    maxEdges(i) = numNeighbors * (numNeighbors - 1);
    combos = combnk(neighborhood, 2);
    neighborhoodPairInds = sub2ind(size(w), [combos(:,1); combos(:,2)], ...
                                            [combos(:,2); combos(:,1)]);
    actualEdges(i) = full(sum(w(neighborhoodPairInds) > 0));                               
  end
  stats.clusteringCoefficient = mean(actualEdges ./ maxEdges);
  
  %
  % All path lengths (full matrix)
  %
  stats.allPathLengths = pathlength(w);
  
  %
  % Characteristic path length
  %
  stats.characteristicPathLength = ...
      full(mean(stats.allPathLengths(~eye(size(w)))));
  
  %
  % Average degree
  %
  stats.averageDegree = full(mean(sum(w'>0)));
  
  %
  % All degrees
  %
  stats.allDegrees = full(sum(w'>0));
  
  %
  % Degree distribution
  %
  stats.degreeDistribution = histc(stats.allDegrees,1:size(w,1));
  
  %
  % PageRanks
  %
  stats.pageRank = pagerank(w');
end
