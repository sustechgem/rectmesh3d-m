% Crop a model/mesh to make it smaller (only keep a trimmed portion)
% FUNCTION [nodeX2, nodeY2, nodeZ2, model2] = cropModel(nodeX1,nodeY1,nodeZ1,model1,xmin,xmax,ymin,ymax,zmax,zmin)
% INPUT
%     nodeX1,nodeY1,nodeZ1,model1: the original mesh and model
%     cropLoc: the volume to be kept, a 6-element vector for the x, y and z ranges
%     [xmin,xmax,ymin,ymax,zmax,zmin]; search for nearest grid; if inf, use extremity
% OUTPUT
%     nodeX2,nodeY2,nodeZ2,model2: the cropped mesh and model
% LAST MODIFIED 20191107 yangdikun@gmail.com
function [nodeX2, nodeY2, nodeZ2, model2] = cropModel(nodeX1,nodeY1,nodeZ1,model1,cropLoc)

Nx = length(nodeX1) - 1;
Ny = length(nodeY1) - 1;
Nz = length(nodeZ1) - 1;

xmin = cropLoc(1);
if isinf(xmin)
    xmin = min(nodeX1);
end

xmax = cropLoc(2);
if isinf(xmax)
    xmax = max(nodeX1);
end

ymin = cropLoc(3);
if isinf(ymin)
    ymin = min(nodeY1);
end

ymax = cropLoc(4);
if isinf(ymax)
    ymax = max(nodeY1);
end

zmax = cropLoc(5);
if isinf(zmax)
    zmax = max(nodeZ1);
end

zmin = cropLoc(6);
if isinf(zmin)
    zmin = min(nodeZ1);
end

[~, xminInd] = min(abs(nodeX1-xmin));
[~, xmaxInd] = min(abs(nodeX1-xmax));
nodeX2 = nodeX1(xminInd:xmaxInd);
[~, yminInd] = min(abs(nodeY1-ymin));
[~, ymaxInd] = min(abs(nodeY1-ymax));
nodeY2 = nodeY1(yminInd:ymaxInd);
[~, zminInd] = min(abs(nodeZ1-zmin));
[~, zmaxInd] = min(abs(nodeZ1-zmax));
nodeZ2 = nodeZ1(zmaxInd:zminInd);

xyz = CellIndex2PointXYZ(nodeX1,nodeY1,nodeZ1,1:Nx*Ny*Nz);
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
ind = x >= nodeX2(1) & x<=nodeX2(end) & y>=nodeY2(1) & y<=nodeY2(end) & z<=nodeZ2(1) & z>=nodeZ2(end);
model2 = model1(ind);

end
