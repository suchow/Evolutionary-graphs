% returns the unique elements of lst that appear at least m times
function n = uniquer(lst,m)
  u = unique(lst);
  for i = 1:length(u)
    n(i) = length(find(lst==u(i))) >= m;
  end
  n = u(n);
end