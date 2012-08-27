% truncated normal distribution
function y = randtn(mu,sd,varargin)
  
  legal = false(varargin{:});
  y = zeros(varargin{:});
  
  pass = false;
  while(~pass)
    illegal = find(~legal);
    y(illegal) = mu + sd*randn(size(illegal));
    legal = y > 0;
    pass = all(legal(:));
  end
end