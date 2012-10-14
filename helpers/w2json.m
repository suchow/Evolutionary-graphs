% takes a graph adjacency matrix and returns the json for a D3 network.
function w2json(w,filename)
  if(nargin < 2)
    filename = 'tmp';
  end
  fid = fopen([filename '.json'],'w');
  numNodes = size(w,1);
  numLinks = sum(sum(w>0));
  
  % print the nodes
  fprintf(fid,'{\n\"nodes\":[\n');
  for i = 1:numNodes
    fprintf(fid,'{\"name\":\"%d\",\"group\":1}',i-1);
    if(i~=numNodes)
      fprintf(fid,',');
    else
      fprintf(fid,'],');
    end
    fprintf(fid,'\n');
  end
  
  % print the links
  [link_i,link_j] = ind2sub(size(w),find(w>0));
  fprintf(fid,'\n');
  fprintf(fid,'\"links\":[\n');
  for i = 1:numLinks
    fprintf(fid, '{\"source\":%d,\"target\":%d,\"value\":1}', link_i(i)-1, link_j(i)-1);
    if(i~=numLinks)
      fprintf(fid,',');
    end
    fprintf(fid,'\n');
  end
  fprintf(fid,']}');
  fclose(fid);
end