% Get volumes of inquired cells
% FUNCTION vol = cellVolume(nodeX,nodeY,nodeZ,cellInd)
% INPUT
%     nodeX,nodeY,nodeZ: a 3D mesh
%     cellInd: index of inquired cell, can be a vector of integer; if empty
%     return all cells
% OUTPUT
%     vol: a vector containing volumns of the inquired cells
% LAST MODIFIED 20191107 yangdikun@gmail.com
function vol = cellVolume(nodeX,nodeY,nodeZ,cellInd)

sizeX = diff(nodeX);
sizeY = diff(nodeY);
sizeZ = -diff(nodeZ);
Ncell = length(sizeX) * length(sizeY) * length(sizeZ);

if isempty(cellInd)
    cellInd = (1:Ncell)';
end

[a, b, c] = meshgrid(sizeX,sizeZ,sizeY);
v = reshape(a,Ncell,[]) .* reshape(b,Ncell,[]) .* reshape(c,Ncell,[]);

vol = abs(v(cellInd));

end




