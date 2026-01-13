% Translate topo in xyz point to topo repesented by mesh model
% FUNCTION topo10 = PointTopo2MeshTopo(topoxyz,nodeX,nodeY,nodeZ)
% INPUT
%     topoxyz: Scattered topography points specified by a 3-column matrix; 
%     each column for X, Y, and Z 
%     nodeX,nodeY,nodeZ: a 3D mesh
% OUTPUT
%     topo10: Topography discretized on the 3D mesh; descrete topography
%     given in the  form of model vector (0 for air and 1 for earth)
% NOTE
%     topoxyz must contain at least three non-collinear points. 
function topo10 = PointTopo2MeshTopo(topoxyz,nodeX,nodeY,nodeZ)


[Nx, Ny, Nz, ~, ~, ~, ~, ~, ~, centerX, centerY, ~] = getMeshPara(nodeX,nodeY,nodeZ);

inMeshInd = PointXYZ2CellIndex(topoxyz,nodeX,nodeY,nodeZ);

if sum(inMeshInd) == 0 % no point in the mesh area
    [xloc, yloc] = meshgrid(centerX,centerY);
    xloc = reshape(xloc',Nx*Ny,1); % tile_center_x vector
    yloc = reshape(yloc',Nx*Ny,1); % tile_center_y vector
    % F = TriScatteredInterp(topoxyz(:,1),topoxyz(:,2),topoxyz(:,3));
    F = scatteredInterpolant(topoxyz(:,1),topoxyz(:,2),topoxyz(:,3));
    tileElevation = F(xloc,yloc);
    tileElevation = reshape(tileElevation,1,[]);
else % at least one point in the mesh area
    topoxyz(inMeshInd==0,:) = [];
    inMeshInd(inMeshInd==0) = [];
    [xyzInd] = GlobalIndex2DirectionalIndex(Nx,Ny,Nz,inMeshInd);
    xInd = xyzInd(:,1);
    yInd = xyzInd(:,2);
    tileInd = (yInd-1) * Nx + xInd;
    temp = sortrows([tileInd topoxyz(:,3)],1);
    tileInd = temp(:,1);
    elevPt = temp(:,2);
    
    newTileIndicator = [true; diff(tileInd)~=0 ];
    NtileIndplus1 = length(tileInd) + 1;
    newTileInd = [find(newTileIndicator); NtileIndplus1 ];
    NnewTileInd = length(newTileInd) - 1;
    tileElevation = zeros(1,Nx*Ny) + NaN;
    
    for p = 1:NnewTileInd
        startInd = newTileInd(p);
        endInd = newTileInd(p+1) - 1;
        tileNumber = tileInd(startInd);
        tileElevation(tileNumber) = median(elevPt(startInd:endInd));
    end
    
    TBDTileInd = isnan(tileElevation);
    if ~isempty(find(TBDTileInd, 1))
        ABDTileInd = ~TBDTileInd;
        [xloc, yloc] = meshgrid(centerX,centerY);
        xloc = reshape(xloc',Nx*Ny,1); % tile_center_x vector
        yloc = reshape(yloc',Nx*Ny,1); % tile_center_y vector
        ABDtileXYElev = [xloc(ABDTileInd) yloc(ABDTileInd) tileElevation(ABDTileInd)'];
        TBDtileXY = [xloc(TBDTileInd) yloc(TBDTileInd)];
        TBDtileElevation = TriScatteredInterpHybrid(ABDtileXYElev,TBDtileXY,'natural','nearest');
        tileElevation(TBDTileInd) = TBDtileElevation;
    end

end

[~, b] = min(abs(  repmat(tileElevation,Nz+1,1)  -  repmat(nodeZ,1,Nx*Ny)  ));
b = b - 1;
Ntile = Nx * Ny;
topo10 = ones(Nz,Ntile);
for p = 1:Ntile
    topo10(1:b(p),p) = 0;
end
topo10 = reshape(topo10,Nx*Ny*Nz,1);

end



% Do TriScatteredInterp using hybrid methods:
% for points inside the convex hull, use linear or natural
% for points outside the convex hull, use nearest (or assign NaN as requested)
function v2 = TriScatteredInterpHybrid(ptm1,ptm2,methodinHull,methodoutHull)
% ptm1: point model 1, 4-column "x y z v" or 3-column "x y v"
% ptm2: locations of point model 2: 3-column "x y z" or 2-column "x y"
% methodinHull: method for ptm2 in hull, 'linear' or 'natural' or 'default'(= [], = 'linear')
% methodoutHull: method for ptm2 outside hull, a scalar to be assigned or 'nearest' or 'default' (= [], = 'nearest')

Ncolptm1 = size(ptm1,2);
if strcmpi(methodinHull,'natural')
    methodinHull = 'natural';
elseif strcmpi(methodinHull,'nearest')
    methodinHull = 'nearest';
else
    methodinHull = 'linear';
end
if ~isscalar(methodoutHull)
    methodoutHull = 'nearest';
end

%v2 = zeros(size(ptm2,1),1);
switch Ncolptm1
    
    case 3 % 2D problem
        x1 = ptm1(:,1);
        y1 = ptm1(:,2);
        v1 = ptm1(:,3);
        x2 = ptm2(:,1);
        y2 = ptm2(:,2);
        F = scatteredInterpolant(x1,y1,v1,methodinHull);
        v2 = F(x2,y2);
        nanInd = isnan(v2);
        if ~isempty(find(nanInd, 1))
            if isscalar(methodoutHull)
                v2(nanInd) = methodoutHull;
            else
                F = scatteredInterpolant([x1; x2(~nanInd)],[y1; y2(~nanInd)],[v1; v2(~nanInd)],'nearest');
                v2(nanInd) = F(x2(nanInd),y2(nanInd));
            end
        end
        
    case 4 % 3D problem
        x1 = ptm1(:,1);
        y1 = ptm1(:,2);
        z1 = ptm1(:,3);
        v1 = ptm1(:,4);
        x2 = ptm2(:,1);
        y2 = ptm2(:,2);
        z2 = ptm2(:,3);
        F = scatteredInterpolant(x1,y1,z1,v1,methodinHull);
        v2 = F(x2,y2,z2);
        nanInd = isnan(v2);
        if ~isempty(find(nanInd, 1))
            if isscalar(methodoutHull)
                v2(nanInd) = methodoutHull;
            else
                F = scatteredInterpolant([x1; x2(~nanInd)],[y1; y2(~nanInd)],[z1; z2(~nanInd)],[v1; v2(~nanInd)],'nearest');
                v2(nanInd) = F(x2(nanInd),y2(nanInd),z2(nanInd));
            end
        end
        
    otherwise
        disp('ERROR: TriScatteredInterpHybrid: wrong # of column of point model 1');
        return;
end


end


