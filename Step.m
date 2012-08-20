% takes one step of a generalized moran process
function population = Step(population)

    % select an individual with probability proportional to fitness
    i = Select(population,0);

    % look at the selected individual's outgoing edges, and choose
    % another individual who this individual's offspring will replace
    j = logical(mnrnd(1,full(population.graph(i,:))));

    % replace the individual
    population.tags(j) = population.tags(i);
    population.memories(j,:) = population.memories(i,:);
end