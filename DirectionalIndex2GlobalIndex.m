% Get global indices of grid points with their directional indices given
% FUNCTION globalInd = DirectionalIndex2GlobalIndex(Nx,Ny,Nz,directionalxyz)
% INPUT
%     Nx, Nz: number of grid points in x and z
%     directionalxyz: 3-column integer matrix, each column is the 
%     directional indices of the grid points [x-index, y-index, z-index]
% OUTPUT
%     globalInd: global indices after reshaping from a 3D array to a vector
% LAST MODIFIED 20210906 yangdikun@gmail.com
function globalInd = DirectionalIndex2GlobalIndex(Nx,Nz,directionalxyz)

globalInd = (directionalxyz(:,2)-1)*Nx*Nz + (directionalxyz(:,1)-1)*Nz + directionalxyz(:,3);

end