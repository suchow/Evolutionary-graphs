function population = Endow(population,K)

    N = size(population.graph,1);

    % pick types (i.e., which of the K events each individual observed)
    population.tags = randi(K,N,1);

    % pick memories. memories are made up of two numbers: an event assignment
    % and a unique identifier.
    population.memories = [population.tags, [1:N]'];
end