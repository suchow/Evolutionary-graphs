% takes in an adjacency
function stats = GraphStatistics(w,nargin)
  stats.C = ClusteringCoefficient(w);
end

% the average local clustering coefficient
function C = ClusteringCoefficient(w)
  for i = 1:length(w) % find the neighbors
    neighborhood = find(w(i,:));
    neighborhood(ismember(neighborhood, i)) = [];
    numNeighbors = length(neighborhood);
    maxEdges(i) = numNeighbors * (numNeighbors - 1);
    combos = combnk(neighborhood, 2);
    neighborhoodPairInds = sub2ind(size(w), [combos(:,1); combos(:,2)], ...
                                            [combos(:,2); combos(:,1)]);
    actualEdges(i) = full(sum(w(neighborhoodPairInds) > 0));                               
  end
  C = mean(actualEdges ./ maxEdges);
end