function [numRemembered, precision] = GetForgettingFunction(history)
  L = 0.74; % from bays and husain
  for i = 1:length(history)
    s = PopulationStatistics(history{i});
    numRemembered(i) = s.numRemembered;
    numSurvivingLineages(i) = s.numSurvivingLineages;
  end
  proportionOfResource = numSurvivingLineages ./ (numRemembered * history{1}.N);
  precision = sqrt(1 ./ (0.004 .* (proportionOfResource .^ L)));
end
