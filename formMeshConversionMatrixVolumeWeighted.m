% Convert model1 defined on mesh1 to model2 defined on mesh2
% Express the conversion as Q * m1 = m2, this function calculate the matrix Q 
% FUNCTION Q = formMeshConversionMatrixVolumeWeighted(nodeX1,nodeY1,nodeZ1,topo1,nodeX2,nodeY2,nodeZ2,topo2,outsideM1method)
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
%     Algorithm based on averaging weighted by overlapping volumes
%     Inactive cells in m2 will be assigned NaN
% LAST MODIFIED 20191112 yangdikun@gmail.com
function Q = formMeshConversionMatrixVolumeWeighted(nodeX1,nodeY1,nodeZ1,topo1,nodeX2,nodeY2,nodeZ2,topo2,outsideM1method)

% mesh prep
Nx1 = length(nodeX1) - 1;
Ny1 = length(nodeY1) - 1;
Nz1 = length(nodeZ1) - 1;
N1 = Nx1*Ny1*Nz1;

Nx2 = length(nodeX2) - 1;
Ny2 = length(nodeY2) - 1;
Nz2 = length(nodeZ2) - 1;
N2 = Nx2*Ny2*Nz2;

% form intersection relationships, containing info about the length of
% intersecting (hx1,hx2), (hy1,hy2), (hz1,hz2)
temp = sortrows([ [nodeX1 ones(Nx1+1,1)]; [nodeX2 zeros(Nx2+1,1)+2]],[1 2]);
dx = diff(temp(:,1));
ind1 = cumsum(temp(:,2)==1);
ind2 = cumsum(temp(:,2)==2);
ind = ind1>=1 & ind1<=Nx1 & ind2>=1 & ind2<=Nx2;
hx1Ind = ind1(ind);
hx2Ind = ind2(ind);
intersechx = dx(ind);
hx1Ind = [ones(hx2Ind(1)-1,1); hx1Ind; zeros(Nx2-hx2Ind(end),1)+Nx1 ];
intersechx = [zeros(hx2Ind(1)-1,1)+999; intersechx; zeros(Nx2-hx2Ind(end),1)+999 ]; % 999 assigned to no-overlapping cells
hx2Ind = [(1:hx2Ind(1)-1)'; hx2Ind; (hx2Ind(end)+1:Nx2)'];
ind = intersechx == 0;
hx1Ind(ind) = [];
hx2Ind(ind) = [];
intersechx(ind) = [];

temp = sortrows([ [nodeY1 ones(Ny1+1,1)]; [nodeY2 zeros(Ny2+1,1)+2]],[1 2]);
dy = diff(temp(:,1));
ind1 = cumsum(temp(:,2)==1);
ind2 = cumsum(temp(:,2)==2);
ind = ind1>=1 & ind1<=Ny1 & ind2>=1 & ind2<=Ny2;
hy1Ind = ind1(ind);
hy2Ind = ind2(ind);
intersechy = dy(ind);
hy1Ind = [ones(hy2Ind(1)-1,1); hy1Ind; zeros(Ny2-hy2Ind(end),1)+Ny1 ];
intersechy = [zeros(hy2Ind(1)-1,1)+999; intersechy; zeros(Ny2-hy2Ind(end),1)+999 ]; % 999 assigned to no-overlapping cells
hy2Ind = [(1:hy2Ind(1)-1)'; hy2Ind; (hy2Ind(end)+1:Ny2)'];
ind = intersechy == 0;
hy1Ind(ind) = [];
hy2Ind(ind) = [];
intersechy(ind) = [];

temp = sortrows([ [nodeZ1 ones(Nz1+1,1)]; [nodeZ2 zeros(Nz2+1,1)+2]],[-1 2]);
dz = -diff(temp(:,1));
ind1 = cumsum(temp(:,2)==1);
ind2 = cumsum(temp(:,2)==2);
ind = ind1>=1 & ind1<=Nz1 & ind2>=1 & ind2<=Nz2;
hz1Ind = ind1(ind);
hz2Ind = ind2(ind);
intersechz = dz(ind);
hz1Ind = [ones(hz2Ind(1)-1,1); hz1Ind; zeros(Nz2-hz2Ind(end),1)+Nz1 ];
intersechz = [zeros(hz2Ind(1)-1,1)+999; intersechz; zeros(Nz2-hz2Ind(end),1)+999 ]; % 999 assigned to no-overlapping cells
hz2Ind = [(1:hz2Ind(1)-1)'; hz2Ind; (hz2Ind(end)+1:Nz2)'];
ind = intersechz == 0;
hz1Ind(ind) = [];
hz2Ind(ind) = [];
intersechz(ind) = [];

% all possible combination of intersecting volumes
[a, b, c] = meshgrid(hx1Ind,hz1Ind,hy1Ind);
hx1hy1hz1Ind = [reshape(a,[],1) reshape(c,[],1) reshape(b,[],1) ];
[a, b, c] = meshgrid(hx2Ind,hz2Ind,hy2Ind);
hx2hy2hz2Ind = [reshape(a,[],1) reshape(c,[],1) reshape(b,[],1) ];
[a, b, c] = meshgrid(intersechx,intersechz,intersechy);
intersechxhyhz = [reshape(a,[],1) reshape(c,[],1) reshape(b,[],1) ];
intersecVolume = prod(intersechxhyhz,2); % get volume = dx * dy * dz

% translate directional cell index to global cell index of mesh
m1Ind = DirectionalIndex2GlobalIndex(Nx1,Ny1,Nz1,hx1hy1hz1Ind);
m2Ind = DirectionalIndex2GlobalIndex(Nx2,Ny2,Nz2,hx2hy2hz2Ind);


%%%% have done universal calculation, next work on special treatment

m1airind = find(topo1==0);
%m1earthind = find(topo1==1);
m2airind = find(topo2==0);
%m2earthind = find(topo2==1);

% remove intersecting volume involving inactive cells in m1 or m2
temp = ismember(m1Ind,m1airind) | ismember(m2Ind,m2airind);
m1Ind(temp) = [];
m2Ind(temp) = [];
intersecVolume(temp) = [];

% remove intersecting volume involving m2 cells outside m1 if outsideM1method = NaN, Inf or 0
m2outInd = false(N2,1);
if ~isempty(outsideM1method)
    directionalInd = GlobalIndex2DirectionalIndex(Nx2,Ny2,Nz2,(1:N2)');
    m2outInd = (nodeX2(directionalInd(:,1)+1) < nodeX1(1) | ...
           nodeX2(directionalInd(:,1)) > nodeX1(end) | ...
           nodeY2(directionalInd(:,2)+1) < nodeY1(1) | ...
           nodeY2(directionalInd(:,2)) > nodeY1(end) | ...
           nodeZ2(directionalInd(:,3)+1) > nodeZ1(1) | ...
           nodeZ2(directionalInd(:,3)) < nodeZ1(end)) & (topo2==1); % active cells in m2 outside of m1 mesh

    temp = ismember(m2Ind,find(m2outInd));
    m1Ind(temp) = [];
    m2Ind(temp) = [];
    intersecVolume(temp) = [];
end 

% normalize weights
sumrow = zeros(N2,1);
Nm2Ind = length(m2Ind);
for p = 1:Nm2Ind
    sumrow(m2Ind(p)) = sumrow(m2Ind(p)) + intersecVolume(p);
end
sumrow(sumrow==0) = 1; % avoid divided by zero
for p = 1:Nm2Ind
    intersecVolume(p) = intersecVolume(p) / sumrow(m2Ind(p));
end

% some active m2 cells (on surface) do not have overlapping with active m1 cells
m2add = logical(topo2) & (~m2outInd);
m2add(m2Ind) = false;
m2addInd = find(m2add)';

addm1Ind = [];
addm2Ind = [];
addintersecVolume = [];
for p = m2addInd
    ind = find(m2Ind>p,1);
    ind = m2Ind==m2Ind(ind);
    addm1Ind = [addm1Ind; m1Ind(ind)];
    addm2Ind = [addm2Ind; zeros(sum(ind),1)+p];
    addintersecVolume = [addintersecVolume; intersecVolume(ind)];
end

% outsider
Nm2outInd = sum(m2outInd);
if isinf(outsideM1method)
    outm1Ind = ones(Nm2outInd,1);
    outm2Ind = find(m2outInd);
    outintersecVolume = zeros(Nm2outInd,1) + Inf;
elseif isnan(outsideM1method)
    outm1Ind = ones(Nm2outInd,1);
    outm2Ind = find(m2outInd);
    outintersecVolume = zeros(Nm2outInd,1) + NaN;
elseif outsideM1method==0
    outm1Ind = ones(Nm2outInd,1);
    outm2Ind = find(m2outInd);
    outintersecVolume = zeros(Nm2outInd,1);
else % empty
    outm1Ind = [];
    outm2Ind = [];
    outintersecVolume = [];
end

% topo2 inactive cells = NaN
airm1Ind = ones(length(m2airind),1);
airm2Ind = m2airind;
airintersecVolume = zeros(length(m2airind),1) + NaN;

% assembly
Q = sparse(m2Ind,m1Ind,intersecVolume,N2,N1) + ...
    sparse(addm2Ind,addm1Ind,addintersecVolume,N2,N1) + ...
    sparse(outm2Ind,outm1Ind,outintersecVolume,N2,N1) + ...
    sparse(airm2Ind,airm1Ind,airintersecVolume,N2,N1);


end









