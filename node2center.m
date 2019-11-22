% Calculate location of cell center along x/y/z direction of mesh
% FUNCTION center = node2center(node)
% INPUT
%     node: a vector of node in a particular direction
% OUTPUT
%     center: a vector of cell centers in the same direction of node
% LAST MODIFIED 20191122 yangdikun@gmail.com
function center = node2center(node)

node = reshape(node,[],1);
Nnode = length(node);
W = spdiags(ones(Nnode-1,2),[0 1],Nnode-1,Nnode);
center = W * node * 0.5;

end