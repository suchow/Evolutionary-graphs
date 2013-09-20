function x = pagerank(w,p)
% PAGERANK  Google's PageRank

if nargin < 3
  p = .85; 
end

w = w - diag(diag(w));
[n,n] = size(w);
c = sum(w,1);
r = sum(w,2);
k = find(c~=0);
D = sparse(k,k,1./c(k),n,n);
e = ones(n,1);
w = p*w*D;
z = ((1-p)*(c~=0) + (c==0))/n;
x = ones(n,1)/n;
xprev = 1;
while sum(abs(xprev-x)) > 0.001
    xprev = x;
    x = w*x + e*(z*x);
end
x = x/sum(x);
