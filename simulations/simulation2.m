%
% Fit thresholded Moran process to data
%
function results = simulation2(nReal)
  Ks = [1,2,3,4,6,8,12];
  numSteps = 2000;
  numIterations = 32;
  ts = [0.0312, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16];
  
  M = 50;
  p = 0.75;
  tau = 0.01;
  threshold = 7;

  results.final = fminsearch(@computeError,[M,p,tau,threshold])
  
  function err = computeError(params)
    params
    N = simulateProcess(params(1),Ks,params(2),numSteps,params(3),numIterations,ts,params(4));
    nReal;
    err = sum(sum((nReal - N).^2))
    
    % plot timecourse of number remembered
    colors = palettablecolors(length(Ks));
    for kIndex = 1:length(Ks)
      subplot(2,1,1)
      plot(ts,N(kIndex,:), '.-', 'Color', colors(kIndex,:));
      hold on;
      axis([0 max(ts) 0 4])
      makepalettable
      
      subplot(2,1,2)
      plot(ts,nReal(kIndex,:), '.-', 'Color', colors(kIndex,:));
      hold on;
      axis([0 max(ts) 0 4])
      makepalettable
    end
    drawnow;
    ylabel('Number of remembered objects');
    xlabel('Time (s)')
    results.nPredicted = N;
    results.nReal = nReal;
    results.SS = sum(sum((nReal - N).^2));
    results.rho = corr([N(:),nReal(:)]);
  end
end


function N = simulateProcess(M,Ks,p,numSteps,tau,numIterations,ts,threshold)
  
  %
  % Run the simulation
  %
  figure(22);
  
  % sanitize input
  M = round(M);
  threshold = round(threshold);
        
  for kIndex = 1:length(Ks);
    K = Ks(kIndex);
    % preallocate for speed
    history = cell(numIterations);
    N = zeros(numSteps+1,numIterations);
    for j = 1:numIterations
      history{j} = ThresholdedMoranProcess2(M,K,p,numSteps,{'Complete'},threshold);
      % compute stats for this iteration
      if(~iscell(history{j}) && (history{j} == -1))
        N(:,j) = zeros(numSteps+1,1);
      else  
        N(:,j) = GetForgettingFunction(history{j});
      end
    end    
    allN(:,kIndex) = mean(N,2);
  end
  totalSimulatedTime = tau*numSteps;
  whichSteps = round((ts./totalSimulatedTime)*numSteps);
  N = allN(whichSteps,:)';
  clf
end
