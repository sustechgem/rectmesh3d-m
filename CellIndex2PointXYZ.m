% Get xyz location of cell centers
% FUNCTION xyz = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,ind)
% INPUT
%     nodeX,nodeY,nodeZ: a 3D mesh
%     ind: cell index number in a mesh, can be a vector; if empty return
%     all cells
% OUTPUT
%     xyz: a 3-column matrix; each row for one cell; columns for x,y,z
% LAST MODIFIED 20191107 yangdikun@gmail.com
function xyz = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,ind)

Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;
centerX = node2center(nodeX);
centerY = node2center(nodeY);
centerZ = node2center(nodeZ);

if isempty(ind)
    ind = (1:Nx*Ny*Nz)';
end

% calculate x y z index
yind = ceil(ind ./ (Nx*Nz));
xind = ceil( (ind-(yind-1).*Nx*Nz) ./ Nz );
zind = ind - (yind-1).*Nx*Nz - (xind-1).*Nz;

% output
xyz = [centerX(xind) centerY(yind) centerZ(zind)];

end

