% removes from lst all elements with fewer than m representatives
function lst = killrare(lst, m)
  unq = unique(lst);
  for i = 1:length(unq)
    matches = (lst==unq(i));
    if(sum(matches) < m)
      lst(matches) = [];
    end
  end
end