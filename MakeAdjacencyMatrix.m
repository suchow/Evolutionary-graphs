% Builds a function that, given indices i,j < N, returns the transition
% probability for a graph of type graphType with N nodes.
function w = MakeAdjacencyMatrix(graphType,N,varargin)
  
  switch graphType
  

  % Complete graph with identically-weighted edges, as in the Moran Process
  case 'Complete'
    w = 1/N .* ones(N);
    
  
  % A complete bipartite graph
  case 'Complete bipartite'
    n = varargin{1};
    m = varargin{2};
    if(n+m ~= N)
      error('In a bipartite graph, the two sides must sum to N.')
    end
    w = sparse(N,N);
    w(1:n,n+1:end) = 1;
    w((n+1):N,1:n) = 1;
    w = rownormalize(w);
  
  % The G(n, p) Erdős–Rényi model proposed in Edgar Gilbert's 1959 "Random
  % graphs", except with no isolated vertices. Defaults to p = 0.5 unless
  % specified as an argument to MakeAdjacencyMatrix.
  case 'Random'
    if(length(varargin) < 1)
      p = 0.5;
    else
      p = varargin{1};
    end
    w = sparse(N,N);
    done = false;
    while(~all(sum(w,2) ~= 0)) % resample all empty rows
      emptyRows = (sum(w,2) == 0);
      w(emptyRows,:) = (rand(sum(emptyRows),N) < p);    
    end  
    w = rownormalize(w);


  % Square lattice connected to its von Neumann neighborhood
  case 'Square lattice 4'
    r = sqrt(N);
    c = sqrt(N);
    diag1 = repmat([ones(c-1,1); 0],r,1);  
    diag1 = diag1(1:end-1);
    diag2 = ones(c*(r-1),1);
    w = diag(diag1,1) + diag(diag2,c);
    w = w+w.';  
    w = rownormalize(w);  


  % Square lattice connected to its Moore neighborhood
  case 'Square lattice 8'
    r = sqrt(N);
    c = sqrt(N);
    diag1 = repmat([ones(c-1,1); 0],r,1);
    diag1 = diag1(1:end-1);            
    diag2 = [0; diag1(1:(c*(r-1)))]; 
    diag3 = ones(c*(r-1),1);         
    diag4 = diag2(2:end-1);       
    w = diag(diag1,1) + diag(diag2,c-1) + ...
        diag(diag3,c) + diag(diag4,c+1);
    w = w+w.';                    
    w = rownormalize(w);  
          

  % A path from vertex 1 to N
  case 'Path'
    w = sparse(N,N);
    w = addlinks(w, 1:N, min(2:(N+1),N));

  % A directed cycle created by having the snake eat its tail
  case 'Cycle'
    w = sparse(N,N);
    w = linktoneighbors(w,1);


  case 'Star'
    w = sparse(N,N);
    w(1,2:end) = true;
    w(2:end,1) = true;
    w = rownormalize(w);


  case 'Burst'
    w = sparse(eye(N));
    w(1,1) = false;
    w(1,2:end) = true;
    w = rownormalize(w);

  % All nodes connect to a central hub
  case 'Integrator'
    w = sparse(N,N);
    w(:,1) = true;


  % Uses the Barabasi-Albert model to construct a scale free network with
  % scaling exponent 3. We seed it with m_0 = 2 nodes, and an edge from 
  % node 1 to node 2 and back. Each step connects a new node to two existing
  % nodes with probability proprtional to their degree
  case 'Scale free'
    w = sparse(N,N);
    w(1,2) = true;
    w(2,1) = true;
    w(3,1) = true;
    w(1,3) = true;
    for i = 4:N
      p_i = rownormalize(sum(full(w)')); % count edges
      n1 = find(mnrnd(1, p_i));          % pick first node
      p_i(n1) = 0;                       % disallow repeat (w/o replacement)
      p_i = rownormalize(p_i);           % renormalize
      n2 = find(mnrnd(1, p_i));          % pick second node

      w(i,n1) = true;                    % link up the chosen nodes
      w(n1,i) = true;
      w(i,n2) = true;
      w(n2,i) = true;
    end
    w = rownormalize(w);
    
  
  case 'Friendship'
    numMills = (N - 1) / 2;
    w = sparse(N,N);
    w(1,2:end) = true;
    for i = 1:numMills
      w(2*i,2*i+1) = true;
      w(2*i+1,2*i) = true;
      w(2*i,1) = true;
      w(1,2*i) = true;
      w(2*i+1,1) = true;
      w(1,2*i+1) = true;
    end
    w = rownormalize(w);
    
  % From Mertzios et al. (2012)  
  case 'Clique wheel'
    w = sparse(N,N);
    n = N/2;
    for i = 1:n % link inner to outer wheel
      w(i,n+i) = true;
      w(n+i,i) = true;
    end
    for i = 1:(n-1)
      w(i,i+1) = true; % inner links
      w(i+1,i) = true;  
      w(n+i,n+i+1) = true; % outer links
      w(n+i+1,n+i) = true;
    end
    w(1,n) = true; % tie together inner
    w(n,1) = true;
    w(N,n+1) = true; % tie together outer
    w(n+1,N) = true;
    w = rownormalize(w);
    
  case 'Konigsberg'
    w = sparse(4,4);
    w(1,2) = 2;
    w(2,1) = 2;
    w(1,3) = 2;
    w(3,1) = 2;
    w(1,4) = 1;
    w(4,1) = 1;
    w(3,4) = 1;
    w(4,3) = 1;
    w(2,4) = 1;
    w(4,2) = 1;
    w = rownormalize(w);

  % Watts and Strogatz model
  case 'Watts-Strogatz'
    if(nargin < 3)
      p = 0.5;
    else
      p = varargin{1}/2;
    end
    w = sparse(N,N);
    w = linktoneighbors(w,[1,2,-1,-2]);
    % random rewiring
    [i,j] = find(w);
    for k = 1:length(i)
      if((rand < p) & (j > i))
        % break the old link between i and j
        w(i(k),j(k)) = false;
        w(j(k),i(k)) = false;
        % choose a new target
        newTarget = randi(N);
        % create link from i to target
        w(i(k),newTarget) = true;
        w(newTarget,i(k)) = true;
      end
    end
    w = rownormalize(w);

  otherwise
    error('This graph type is not supported.')
    
  end
end

function y = rownormalize(x)
  y = bsxfun(@times, x, 1./(sum(x, 2)));
end

function w = addlinks(w,from,to)
  w(sub2ind(size(w),from,to)) = true;
end

function w = linktoneighbors(w,t)
  N = length(w);
  w = addlinks(w,1:N,1+mod(t+[0:(N-1)],N));
  for i = 1:length(t)
    w = addlinks(w,1:N,1+mod(t(i)+[0:(N-1)],N));
  end
end