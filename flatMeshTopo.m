% Make a flat topography model
% FUNCTION topo = flatMeshTopo(nodeX,nodeY,nodeZ,elev)
% INPUT
%     nodeX,nodeY,nodeZ: a 3D mesh
%     elev: a scalar specifying the elevation
% OUTPUT
%     topo: a topography model vector; 0 = above topo; 1 = below topo
% LAST MODIFIED 20191107 yangdikun@gmail.com
function topo = flatMeshTopo(nodeX,nodeY,nodeZ,elev)

[~, b] = min(abs(nodeZ-elev));
b = b - 1;

onecolumn = ones(length(nodeZ)-1,1);
onecolumn(1:b) = 0;

topo = repmat(onecolumn, (length(nodeX)-1)*(length(nodeY)-1) ,1);


end