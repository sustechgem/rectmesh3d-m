% Get directional index of cells in a mesh with global cell index given
% FUNCTION directionalInd = GlobalIndex2DirectionalIndex(Nx,Ny,Nz,globalInd)
% INPUT
%     Nx, Ny, Nz: number of cells in x, y, z
%     globalInd: global cell index in the model vector
% OUTPUT
%     directionalxyz: 3-column integer matrix, each column is for the 
%     directional index of the inquiry cells (1st column: index counting 
%     in +x direction; 2nd column: index counting in +y direction; 3rd
%     column: index counting in -z direction.
% NOTE
%     Because of how the cells are counted globally, Ny is not actually
%     useful, but included as input parameter for completeness
% LAST MODIFIED 20191122 yangdikun@gmail.com
function directionalInd = GlobalIndex2DirectionalIndex(Nx,Ny,Nz,globalInd)

globalInd = reshape(globalInd,[],1);
yind = ceil(globalInd./Nx./Nz);
xind = ceil( (globalInd-(yind-1).*Nx.*Nz) ./ Nz );
zind = globalInd - (yind-1).*Nx.*Nz - (xind-1).*Nz;
directionalInd = [xind yind zind];
end



