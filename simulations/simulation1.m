%
% Simuates a neutral memory process on various graph structures, plotting
% the forgetting function and fidelity function at various set sizes.
%
function results = simulation1

  N = 63;
  Ks = [1,2,3,4];
  numSteps = 1000;
  numIterations = 8;

  graphTypes = {'Complete'};
  % graphTypes = {'Complete', {'Visual',N,1}, 'Square lattice 4', 'Cycle', ...
  %                'Scale free', 'Erdos-Renyi', 'Square lattice 8', 'Star', ...
  %                'Burst', 'Integrator', 'Model A', 'Clique wheel', ...
  %                'Bethe lattice', 'Random regular'}

  %
  % Run the simulation
  %
  figure;
  for q = 1:length(graphTypes)
    if(~iscell(graphTypes{q}))
      thisGraphType = {graphTypes{q}};
    else
      thisGraphType = graphTypes{q};
    end

    for kIndex = 1:length(Ks);
      K = Ks(kIndex);
      % preallocate for speed
      history = cell(numIterations);
      numRemembered = zeros(numSteps+1,numIterations);
      precision = zeros(numSteps+1,numIterations);
      for j = 1:numIterations
        history{j} = ThresholdedMoranProcess(N,K,numSteps,thisGraphType,16);
        % compute stats for this iteration
        [numRemembered(:,j), precision(:,j)] = GetForgettingFunction(history{j});
      end

      results.mean(q) = mean(numRemembered(end,:));
      results.sd(q) = std(numRemembered(end,:))./sqrt(numIterations)

      fprintf('%s\n', thisGraphType{1})
      fprintf('events: %2.2f, precision: %2.2f\n\n', ...
              mean(numRemembered(end,:)), ...
              mean(precision(end,:)))

      colors = palettablecolors(length(Ks));

      % plot timecourse of number remembered
      subplot(2,length(graphTypes),q)
      plot(mean(numRemembered,2), 'Color', colors(kIndex,:), ...
          'LineWidth', 1.5, 'LineSmoothing', 'off');
      hold on;
      axis([1 numSteps 0 max(Ks)])
      title(graphTypes{q}, 'FontSize', 10);
      if(q == 1), ylabel('Number of remembered events'); end
      makepalettable

      % plot timecourse of
      subplot(2,length(graphTypes),length(graphTypes)+q)
      plot(mean(precision,2), 'Color', colors(kIndex,:), ...
          'LineWidth', 1.5, 'LineSmoothing', 'off');
      hold on;
      drawnow
      makepalettable
      if(q == 1), ylabel('Number of unique observations'); end
      axis([1 numSteps 0 100])
      xlabel('Number of steps')
    end
  end
end
