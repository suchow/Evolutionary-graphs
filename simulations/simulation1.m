%
% Simuates a neutral memory process on various graph structures.
%
function simulation1

  % population settings
  N = 64;
  Ks = [1,2,3,6,12];
  numSteps = 2000;
  graphTypes = {'Complete', 'Square lattice 4', 'Square lattice 8', 'Scale free'} 
  %graphTypes = {'Complete', 'Random', 'Square lattice 4', 'Path', ...
  %    'Cycle', 'Scale free', 'Burst', 'Star', 'Integrator', ...
  %    'Watts-Strogatz', 'Konigsberg', 'Clique wheel'};

  %graphTypes = {'Complete', 'Random', 'Scale free'}
  %graphTypes = {'Bethe lattice', 'Scale free'}

  % simulation settings
  numIterations = 32;

  %
  % Run the simulation
  %
  for q = 1:length(graphTypes)
    for kIndex = 1:length(Ks);
      K = Ks(kIndex);
      numRemembered = zeros(numSteps,numIterations);
      numSurvivingLineages = zeros(numSteps,numIterations);
      for j = 1:numIterations
          population = MakePopulation(N,graphTypes{q});
          population = Endow(population,K);

          for i = 1:numSteps+1
              numRemembered(i,j) = length(unique(population.tags));
              numSurvivingLineages(i,j) = length(unique(population.memories(:,2)));
              population = Step(population);
          end
      end
        
      fprintf('%s\n', graphTypes{q})
      fprintf('events: %2.2f, lineages: %2.2f\n\n', mean(numRemembered(end,:)), ...
                        mean(numSurvivingLineages(end,:)))                          
                  
      colors = palettablecolors(length(Ks));

      % plot timecourse of number remembered
      subplot(2,length(graphTypes),q)
      plot(mean(numRemembered'), 'Color', colors(kIndex,:), ...
          'LineWidth', 1.5, 'LineSmoothing', 'off');
      hold on;
      axis([1 numSteps 0 max(Ks)])
      title(graphTypes{q}, 'FontSize', 16);
      if(q == 1) ylabel('Number of remembered events'); end
      makepalettable

      % plot timecourse of 
      subplot(2,length(graphTypes),length(graphTypes)+q)
      plot(mean(numSurvivingLineages'), 'Color', colors(kIndex,:), ...
          'LineWidth', 1.5, 'LineSmoothing', 'off');
      hold on;
      drawnow
      makepalettable
      if(q == 1) ylabel('Number of unique observations'); end
      axis([1 numSteps 0 N])
      xlabel('Number of steps')
    end
  end
end

% function y = unique(x)
%     x = sort(x);
%     y = x([true;diff(x(:))>0]);
% end
