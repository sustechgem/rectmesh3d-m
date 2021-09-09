% Form an interpolation matrix that interpolates values on a lattice grid to
% arbitrary points using trilinear interpolation
% FUNCTION P = formLatticeTrilinearInterpMatrix(gridX,gridY,gridZ,points)
% INPUT
%     gridX,gridY,gridZ: define a lattice point structure (eg. nodes of a 3D mesh)
%     points: 3-column matrix for the x-y-z locations of the inquiry points
% OUTPUT
%     P: a sparse projection matrix; this matrix can be applied to a vector
%     holding values defined on the lattice grid points
% NOTE
%     gridX,gridY,gridZ can define any lattice structure, for example, can 
%     be the midpoints of edges or the centers of faces or nodes in a 3D
%     mesh.
% LAST MODIFIED 20210906 yangdikun@gmail.com
function P = formLatticeTrilinearInterpMatrix(gridX,gridY,gridZ,points)

Nnode = length(gridX) * length(gridY) * length(gridZ);
Np = size(points,1);
x = points(:,1);
y = points(:,2);
z = points(:,3);
Nx = length(gridX)-1;
Ny = length(gridY)-1;
Nz = length(gridZ)-1;

% in case a point is beyond lattice limits, find the nearest for them
% to make sure all inquiry points are within the lattice structure
x(x<gridX(1)) = gridX(1);
x(x>gridX(end)) = gridX(end);
y(y<gridY(1)) = gridY(1);
y(y>gridY(end)) = gridY(end);
z(z>gridZ(1)) = gridZ(1);
z(z<gridZ(end)) = gridZ(end);

% Trilinear interp: a point in a cubic volume; the weight of a particular
% vertex is proportional to its cooresponding 3D diagonally opposite volume.

% convert points to cellInd, then to directional ind
cellInd = PointXYZ2CellIndex([x y z],gridX,gridY,gridZ);
directionalInd = GlobalIndex2DirectionalIndex(Nx,Ny,Nz,cellInd);
xind = directionalInd(:,1);
yind = directionalInd(:,2);
zind = directionalInd(:,3); 
% nodes ind for the enclosing cube (used for entry position in projection matrix)
n1ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind yind zind]); % node # 1
n2ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind yind zind+1]); % node # 2
n3ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind+1 yind zind]); % node # 3
n4ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind+1 yind zind+1]); % node # 4
n5ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind yind+1 zind]); % node # 5
n6ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind yind+1 zind+1]); % node # 6
n7ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind+1 yind+1 zind]); % node # 7
n8ind = DirectionalIndex2GlobalIndex(Nx+1,Ny+1,Nz+1,[xind+1 yind+1 zind+1]); % node # 8
% location of the enclosing cube (used for weights in projection matrix)
xmin = gridX(xind); xmax = gridX(xind+1);
ymin = gridY(yind); ymax = gridY(yind+1);
zmax = gridZ(zind); zmin = gridZ(zind+1);
dx1 = x-xmin; dx2 = xmax-x; % sub-cell dimensions
dy1 = y-ymin; dy2 = ymax-y; % sub-cell dimensions
dz1 = zmax-z; dz2 = z-zmin; % sub-cell dimensions
vol = (xmax-xmin) .* (ymax-ymin) .* (zmax-zmin); % cell volumes
n8wgt = dx1 .* dy1 .* dz1 ./ vol; % normalized vol of sub-cell # 1 = weight for node # 8
n7wgt = dx1 .* dy1 .* dz2 ./ vol; % normalized vol of sub-cell # 2 = weight for node # 7
n6wgt = dx2 .* dy1 .* dz1 ./ vol; % normalized vol of sub-cell # 3 = weight for node # 6
n5wgt = dx2 .* dy1 .* dz2 ./ vol; % normalized vol of sub-cell # 4 = weight for node # 5
n4wgt = dx1 .* dy2 .* dz1 ./ vol; % normalized vol of sub-cell # 5 = weight for node # 4
n3wgt = dx1 .* dy2 .* dz2 ./ vol; % normalized vol of sub-cell # 6 = weight for node # 3
n2wgt = dx2 .* dy2 .* dz1 ./ vol; % normalized vol of sub-cell # 7 = weight for node # 2
n1wgt = dx2 .* dy2 .* dz2 ./ vol; % normalized vol of sub-cell # 8 = weight for node # 1
% form projection matrix
I = repmat((1:Np)',8,1);
J = [n1ind; n2ind; n3ind; n4ind; n5ind; n6ind; n7ind; n8ind];
S = [n1wgt; n2wgt; n3wgt; n4wgt; n5wgt; n6wgt; n7wgt; n8wgt]; 
P = sparse(I,J,S,Np,Nnode); % trilinear interp.

end