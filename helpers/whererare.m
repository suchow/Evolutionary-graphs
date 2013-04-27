function loc = whererare(lst, m)
	loc = zeros(size(lst));
	unq = unique(lst);
	for i = 1:length(unq)
	  matches = (lst==unq(i));
    % size(matches)
    % size(loc)
    % '--'
	  if(sum(matches) < m)
	    loc = loc + matches;
    end
  end
  loc = logical(loc);
end