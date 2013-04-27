% count how many unique elements of lst appear at least m times
function n = uniquem(lst,m)
  u = unique(lst);
  for i = 1:length(u)
    n(i) = length(find(lst==u(i))) >= m;
  end
  n = sum(n);
end