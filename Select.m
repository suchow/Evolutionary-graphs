% frequency dependent selection with parameter tau that determines the 
% strength of selection. tau is not yet implemented.
function i = Select(population, tau)

    if(nargin < 2)
        tau = 0;
    end

    if(tau ~= 0)
        warning(['Frequency-dependent selection is not yet availble. ' ...
                 'Defaulting to neutral selection.'])
        tau = 0;
    end

    if(tau == 0) % neutral selection
        i = randi(size(population.graph,1));
    end
end