% Calculate all the parameters of a rectilinear mesh from node info
% FUNCTION [Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, centerX, centerY, centerZ] = getMeshPara(nodeX,nodeY,nodeZ)
% INPUT
%     nodeX,nodeY,nodeZ: a 3D mesh
% OUTPUT
%     Nx, Ny, Nz: number of nodes in x, y, z
%     x0, y0, z0: the coordinate of the first counting node in x, y, z
%     hx, hy, hz: a vector of cell sizes in x, y, z
%     centerX, centerY, centerZ: a vector of coordinates of the cell centers in x, y, z
% LAST MODIFIED 20191122 yangdikun@gmail.com
function [Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, centerX, centerY, centerZ] ...
    = getMeshPara(nodeX,nodeY,nodeZ)

Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

x0 = nodeX(1);
y0 = nodeY(1);
z0 = nodeZ(1);

hx = node2size(nodeX);
hy = node2size(nodeY);
hz = node2size(nodeZ);

centerX = node2center(nodeX);
centerY = node2center(nodeY);
centerZ = node2center(nodeZ);

end