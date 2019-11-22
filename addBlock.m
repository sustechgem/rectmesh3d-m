% Add blocks to model
% FUNCTION model_out = addBlock(nodeX,nodeY,nodeZ,model_in,blkLocInfo,blkValInfo)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter
%     model_in: input model, to which new blocks are added, if omitted, assign
%     0 everywhere in the mesh; can be a scalar
%     blkLocInfo: block location info, Nblk x 6 matrix, for each row [xmin xmax ymin ymax zmax zmin]
%     blkValInfo: block's value, Nblk x 1 vector, each element for one block
% OUTPUT
%     model_out: model vector with objects embedded
% LAST MODIFIED 20191107 yangdikun@gmail.com
function model_out = addBlock(nodeX,nodeY,nodeZ,model_in,blkLocInfo,blkValInfo)

Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

Nblk = size(blkLocInfo,1);
blkLocInfo(blkLocInfo(:,1)==-inf,1) = nodeX(1);
blkLocInfo(blkLocInfo(:,2)==inf,2) = nodeX(end);
blkLocInfo(blkLocInfo(:,3)==-inf,3) = nodeY(1);
blkLocInfo(blkLocInfo(:,4)==inf,4) = nodeY(end);
blkLocInfo(blkLocInfo(:,5)==inf,5) = nodeZ(1);
blkLocInfo(blkLocInfo(:,6)==-inf,6) = nodeZ(end);

for p = 1:Nblk

    xmin = blkLocInfo(p,1);
    xmax = blkLocInfo(p,2);
    ymin = blkLocInfo(p,3);
    ymax = blkLocInfo(p,4);
    zmax = blkLocInfo(p,5);
    zmin = blkLocInfo(p,6);
    
    [~, xminInd] = min(abs(nodeX-xmin));
    [~, xmaxInd] = min(abs(nodeX-xmax));
    nodeX2 = nodeX(xminInd:xmaxInd);
    [~, yminInd] = min(abs(nodeY-ymin));
    [~, ymaxInd] = min(abs(nodeY-ymax));
    nodeY2 = nodeY(yminInd:ymaxInd);
    [~, zminInd] = min(abs(nodeZ-zmin));
    [~, zmaxInd] = min(abs(nodeZ-zmax));
    nodeZ2 = nodeZ(zmaxInd:zminInd);
    
    xyz = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,1:Nx*Ny*Nz);
    x = xyz(:,1);
    y = xyz(:,2);
    z = xyz(:,3);
    ind = x >= nodeX2(1) & x<=nodeX2(end) & y>=nodeY2(1) & y<=nodeY2(end) & z<=nodeZ2(1) & z>=nodeZ2(end);
    model_out(ind) = blkValInfo(p);
    
end


end






