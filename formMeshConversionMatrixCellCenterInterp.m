% Convert model1 defined on mesh1 to model2 defined on mesh2
% Express the conversion as Q * m1 = m2, this function calculate the matrix Q 
% FUNCTION Q = formMeshConversionMatrixCellCenterInterp(nodeX1,nodeY1,nodeZ1,topo1,nodeX2,nodeY2,nodeZ2,topo2,outsideM1method)
% INPUT
%     nodeX1,nodeY1,nodeZ1: define mesh1
%     topo1: a model vector of 1/0 indicating active/inactive cells in mesh1; often
%     used for topography info
%     nodeX2,nodeY2,nodeZ2: define mesh2
%     topo2: a model vector of 1/0 indicating active/inactive cells in mesh2; often
%     used for topography info
%     outsideM1method: options if some cells in m2 cells are outside of m1 mesh. Can be 0,
%     Inf or NaN. If empty([]), use nearest extrapolation
% OUTPUT
%     Q: the mapping matrix so Q * m1 = m2
% NOTE
%     Algorithm based on interpolation of cell centers
%     Inactive cells in m2 will be assigned NaN
% LAST MODIFIED 20191121 yangdikun@gmail.com
function Q = formMeshConversionMatrixCellCenterInterp(nodeX1,nodeY1,nodeZ1,topo1,nodeX2,nodeY2,nodeZ2,topo2,outsideM1method)
% nodeX1,nodeY1,nodeZ1: mesh parameters for model m1
% topo1: mnz topo info for mesh1
% nodeX1,nodeY1,nodeZ1: mesh parameters for model m2
% topo2: mnz topo info for mesh2
% outsideM1method: options for the m2 cells outside of m1 mesh. Can be 0,
%                  Inf or NaN. If empty([]), means extend out (nearest extrapolation).

% mesh prep
Nx1 = length(nodeX1) - 1;
Ny1 = length(nodeY1) - 1;
Nz1 = length(nodeZ1) - 1;
N1 = Nx1*Ny1*Nz1;
centerX1 = node2center(nodeX1);
centerY1 = node2center(nodeY1);
centerZ1 = node2center(nodeZ1);

Nx2 = length(nodeX2) - 1;
Ny2 = length(nodeY2) - 1;
Nz2 = length(nodeZ2) - 1;
N2 = Nx2*Ny2*Nz2;
centerX2 = node2center(nodeX2);
centerY2 = node2center(nodeY2);
centerZ2 = node2center(nodeZ2);

% get distance info in x, y, z direction
dx = bsxfun(@minus,centerX1',centerX2);
indxw = sum(bsxfun(@le,dx,0),2);
indxe = indxw + 1; % index of eastern neighbor m1-point
temp = indxw==0; % some do not have a western neighbor, d set to 999
indxw(temp) = 1;
dxw = abs(dx((indxw-1)*Nx2+(1:Nx2)'));
dxw(temp) = 999;
temp = indxe>Nx1; % some do not have an eastern neighbor, d set to 999
indxe(temp) = Nx1;
dxe = abs(dx((indxe-1)*Nx2+(1:Nx2)'));
dxe(temp) = 999;

dy = bsxfun(@minus,centerY1',centerY2); % distance of each point in m2 from every point in m1 (y direction), each column for each point in m2
indys = sum(bsxfun(@le,dy,0),2); % index of southern neighbor m1-point of each m2-point; zero means no southern neighbor
indyn = indys + 1; % index of norhtern neighbor m1-point
temp = indys==0; % some do not have a southern neighbor, d set to 999
indys(temp) = 1;
dys = abs(dy((indys-1)*Ny2+(1:Ny2)'));
dys(temp) = 999;
temp = indyn>Ny1; % some do not have a northern neighbor, d set to 999
indyn(temp) = Ny1;
dyn = abs(dy((indyn-1)*Ny2+(1:Ny2)'));
dyn(temp) = 999;

dz = bsxfun(@minus,centerZ1',centerZ2); % distance of each point in m2 from every point in m1 (z direction), each column for each point in m2
indzu = sum(bsxfun(@ge,dz,0),2); % index of upper neighbor m1-point of each m2-point; zero means no upper neighbor
indzb = indzu + 1; % index of bottom neighbor m1-point
temp = indzu==0; % some do not have an upper neighbor, d set to 999
indzu(temp) = 1;
dzu = abs(dz((indzu-1)*Nz2+(1:Nz2)'));
dzu(temp) = 999;
temp = indzb>Nz1; % some do not have a bottom neighbor, d set to 999
indzb(temp) = Nz1;
dzb = abs(dz((indzb-1)*Nz2+(1:Nz2)'));
dzb(temp) = 999;

% get m2-point's 8 corner point index in m1-x, y, z direction 
xyzIndCorner1 = [ repmat( reshape(repmat(indxw,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indys,1,Nz2*Nx2)',[],1) ...
                  repmat(indzu,Nx2*Ny2,1) ];
xyzIndCorner2 = [ repmat( reshape(repmat(indxw,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indys,1,Nz2*Nx2)',[],1) ...
                  repmat(indzb,Nx2*Ny2,1) ];
xyzIndCorner3 = [ repmat( reshape(repmat(indxe,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indys,1,Nz2*Nx2)',[],1) ...
                  repmat(indzu,Nx2*Ny2,1) ];              
xyzIndCorner4 = [ repmat( reshape(repmat(indxe,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indys,1,Nz2*Nx2)',[],1) ...
                  repmat(indzb,Nx2*Ny2,1) ];              
xyzIndCorner5 = [ repmat( reshape(repmat(indxw,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indyn,1,Nz2*Nx2)',[],1) ...
                  repmat(indzu,Nx2*Ny2,1) ];              
xyzIndCorner6 = [ repmat( reshape(repmat(indxw,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indyn,1,Nz2*Nx2)',[],1) ...
                  repmat(indzb,Nx2*Ny2,1) ];              
xyzIndCorner7 = [ repmat( reshape(repmat(indxe,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indyn,1,Nz2*Nx2)',[],1) ...
                  repmat(indzu,Nx2*Ny2,1) ];              
xyzIndCorner8 = [ repmat( reshape(repmat(indxe,1,Nz2)',[],1) , Ny2 ,1 ) ...
                  reshape(repmat(indyn,1,Nz2*Nx2)',[],1) ...
                  repmat(indzb,Nx2*Ny2,1) ];              

% corner point's global index in m1
globalIndCorner = zeros(N2,8);
globalIndCorner(:,1) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner1);
globalIndCorner(:,2) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner2);
globalIndCorner(:,3) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner3);
globalIndCorner(:,4) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner4);
globalIndCorner(:,5) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner5);
globalIndCorner(:,6) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner6);
globalIndCorner(:,7) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner7);
globalIndCorner(:,8) = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,xyzIndCorner8);

% calculate un-normalized weight
weightCorner = zeros(N2,8);
weightCorner(:,1) = prod([ repmat(reshape(repmat(dxe,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dyn,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzb,Nx2*Ny2,1)] ,2);
weightCorner(:,2) = prod([ repmat(reshape(repmat(dxe,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dyn,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzu,Nx2*Ny2,1)] ,2);
weightCorner(:,3) = prod([ repmat(reshape(repmat(dxw,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dyn,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzb,Nx2*Ny2,1)] ,2);
weightCorner(:,4) = prod([ repmat(reshape(repmat(dxw,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dyn,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzu,Nx2*Ny2,1)] ,2);
weightCorner(:,5) = prod([ repmat(reshape(repmat(dxe,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dys,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzb,Nx2*Ny2,1)] ,2);
weightCorner(:,6) = prod([ repmat(reshape(repmat(dxe,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dys,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzu,Nx2*Ny2,1)] ,2);
weightCorner(:,7) = prod([ repmat(reshape(repmat(dxw,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dys,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzb,Nx2*Ny2,1)] ,2);
weightCorner(:,8) = prod([ repmat(reshape(repmat(dxw,1,Nz2)',[],1),Ny2,1) ...
                           reshape(repmat(dys,1,Nz2*Nx2)',[],1)  ...
                           repmat(dzu,Nx2*Ny2,1)] ,2);
                   
% deal with topo
topo1airind = find(topo1==0);
temp = ismember(globalIndCorner,topo1airind);
weightCorner(temp) = 0; % 0-weight on inactive cells of m1

% discretized topo of m1 and m2 could be a couple of cell different
% there may be some active cells in m2 that request m1 cells in air
% those cells are re-assigned the same value as surface cells (search down)
ind = flipud(find(sum(temp,2)==8 & topo2==1))'; % m2 cells that have no contributing m1 cells and is active; flipud for backward tracking (bottom-up)
for p = ind
    globalIndCorner(p,:) = globalIndCorner(p+1,:); % copy from the cell below
    weightCorner(p,:) = weightCorner(p+1,:); % copy from the cell below
end 

% normalize weights
sumrow = sum(weightCorner,2);     
sumrow(sumrow==0) = 1; % avoid divided by zero
weightCorner = bsxfun(@rdivide,weightCorner,sumrow); % weightCorner(:,?) ./ sumrow

% if m2 cells outside of m1 mesh
if ~isempty(outsideM1method)
    
    m2xyz = CellIndex2PointXYZ(nodeX2,nodeY2,nodeZ2,(1:N2)');
    outInd = m2xyz(:,1)<centerX1(1) | m2xyz(:,1)>centerX1(end) | ...
             m2xyz(:,2)<centerY1(1) | m2xyz(:,2)>centerY1(end) | ...
             m2xyz(:,3)>centerZ1(1) | m2xyz(:,3)<centerZ1(end);

    if isinf(outsideM1method)
        weightCorner(outInd,1) = Inf;
    elseif isnan(outsideM1method)
        weightCorner(outInd,1) = NaN;
    else % outsideM1method = 0
        weightCorner(outInd,:) = 0;
    end
    
end

weightCorner(topo2==0,1) = NaN; % NaN for inactive cell in m2

% assemble Q matrix
Q = spalloc(N2,N1,N2*8);
for p = 1:8
    Q = Q + sparse((1:N2),globalIndCorner(:,p),weightCorner(:,p),N2,N1,N2*8);
end



end

