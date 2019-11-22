% Get global index of cells in a mesh with directional index given
% FUNCTION globalInd = DirectionalIndex2GlobalIndex(Nx,Ny,Nz,directionalxyz)
% INPUT
%     Nx, Ny, Nz: number of cells in x, y, z
%     directionalxyz: 3-column integer matrix, each column is the 
%     directional index of the cell (i-th row, j-th column, k-th layer)
% OUTPUT
%     globalInd: global cell index in the model vector
% NOTE
%     Because of how the cells are counted globally, Ny is not actually
%     useful, but included as input parameter for completeness
% LAST MODIFIED 20191107 yangdikun@gmail.com
function globalInd = DirectionalIndex2GlobalIndex(Nx,Ny,Nz,directionalxyz)

globalInd = (directionalxyz(:,2)-1)*Nx*Nz + (directionalxyz(:,1)-1)*Nz + directionalxyz(:,3);

end