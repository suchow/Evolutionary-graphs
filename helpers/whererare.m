% Returns a logical array whose ith element is true iff the ith element of lst
% is represented fewer than m times. In other words, it finds the rare elements
% lst, with m as the criterion definition of rarity.
function loc = whererare(lst, m)
	loc = zeros(size(lst));
	unq = unique(lst);
	for i = 1:length(unq)
	  matches = (lst==unq(i));
	  if(sum(matches) < m)
	    loc = loc + matches;
    end
  end
  loc = logical(loc);
end
