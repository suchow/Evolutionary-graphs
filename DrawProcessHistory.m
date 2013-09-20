function DrawProcessHistory(history)
  numSteps = length(history)-1;
  allTypes = unique(history{1}.tags);
  numTypes = length(allTypes);
  colors = palettablecolors(numTypes);
  initialN = history{1}.N;
  xscale = 0.14;
  yscale = 0.15;
  markerSize = 20;
  fig = figure('Visible', 'off', 'PaperSize', [initialN*xscale,numSteps*yscale]);
  for i = 1:(numSteps+1)
    N = history{i}.N;
    % draw the population
    for j = 1:N
      thisTag = history{i}.tags(j);
      if(thisTag == -1)
        color = [0.8,0.8,0.8];
      else
        thisTypeNumber = find(allTypes == history{i}.tags(j));
        color = colors(thisTypeNumber,:);
      end
      plot(j, -i, '.', 'MarkerSize',markerSize,'Color',color);
      hold on;
    end
    % draw all deaths
    death = history{i}.death;
    for j = 1:length(death)
      plot(death(j), -i, '.', 'MarkerSize',markerSize-(markerSize*0.5),'Color',[1,1,1]);
    end
    % draw all reproductions
    reproduce = history{i}.reproduce;
    for j = 1:length(reproduce)
      plot(reproduce(j), -i, '+', 'MarkerSize',2.5,'Color',[1,1,1]);
    end
  end
  grid off
  
  % pretty things up
  makepalettable
  set(gca,'Visible','off')
  
  % output
  set(gcf, 'PaperPositionMode', 'manual');
  set(fig,'PaperPosition',[0 0 initialN*xscale numSteps*yscale])
  print(fig,'-dpdf','test')
  
end