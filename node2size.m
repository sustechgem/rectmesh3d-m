% Calculate size of cell along x/y/z direction of mesh
% FUNCTION h = node2size(node)
% INPUT
%     node: a vector of node in a particular direction
% OUTPUT
%     h: a vector of cell size in the same direction of node
% LAST MODIFIED 20191122 yangdikun@gmail.com
function h = node2size(node)

h = abs(diff(node));

end






