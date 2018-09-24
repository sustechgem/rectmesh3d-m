% Visualize 3D model on a regular mesh. Can display cell values on inquiry
% blocks or sections/slices
% FUNCTION hp = viewModel(nodeX,nodeY,nodeZ,model,logflag,block)
% INPUT
%     nodeX, nodeY, nodeZ: define a 3D mesh
%     model: a vector of model in UBC-GIF ordering; if empty, only load
%         mesh; use NaN to mask the ignored value.
%     logflag: empty or 'log', indicating the color map is in linear or logarithmic scale
%     block: [Nblk x 6] numeric array, specify the boundaries of 
%         blocks to be viewed; each row for a block and can have multiple
%         blocks; within each row [xmin xmax ymin ymax zmax zmin]; use -inf
%         or inf to indicate range of infinity; block boundaries in meter are 
%         snapped to the nearest nodes; if max and min are snapped to the same
%         node, use the mean of max and min relative to the node to determine which side of 
%         the face is taken for the cross section; empty means showing the entire.
% OUTPUT
%    hp: [Nblk x 1] Patch array for the handles of the patch objects
%    created in the axis
function hp = viewModelBlock(nodeX,nodeY,nodeZ,model,logflag,block)



Nnx = length(nodeX);
Nny = length(nodeY);
Nnz = length(nodeZ);
Ncx = Nnx-1;
Ncy = Nny-1;
Ncz = Nnz-1;

         
          
if isempty(model)
    model = nan(Ncx*Ncy*Ncz,1);
end

if strcmpi(logflag,'log')
    model = log10(model);
end

if isempty(block)
    block = [nodeX(1) nodeX(end) nodeY(1) nodeY(end) nodeZ(1) nodeZ(end)];
else
    % deal with inf
    block(block(:,1)==-inf,1) = nodeX(1);
    block(block(:,2)==inf,2) = nodeX(end);
    block(block(:,3)==-inf,3) = nodeY(1);
    block(block(:,4)==inf,4) = nodeY(end);
    block(block(:,5)==inf,5) = nodeZ(1);
    block(block(:,6)==-inf,6) = nodeZ(end);
end
Nblk = size(block,1);

% get blocks' walls
Wall = cell(Nblk,1);
WallView = cell(Nblk,1);
for p = 1:Nblk
    
    % snap to nearest node and convert meter to node count
    [~, wNodeInd] = min(abs(nodeX-block(p,1)));
    [~, eNodeInd] = min(abs(nodeX-block(p,2)));
    [~, sNodeInd] = min(abs(nodeY-block(p,3)));
    [~, nNodeInd] = min(abs(nodeY-block(p,4)));
    [~, tNodeInd] = min(abs(nodeZ-block(p,5)));
    [~, bNodeInd] = min(abs(nodeZ-block(p,6)));
    
    % west and east wall
    [a, b, c]=meshgrid([wNodeInd eNodeInd],sNodeInd:nNodeInd-1,tNodeInd:bNodeInd-1);
    weWall = [a(:) b(:) c(:)];
    if wNodeInd == eNodeInd % for cross section, sign determined by the inquiry location relative to the node
        weWallView = a(:)*0 + sign(mean(block(p,1:2)) - nodeX(wNodeInd) - 1e-10); % in case inquiry location exactly on the node, give a slight negative push so east wall dominates
    else
        weWallView = a(:)*0;
        weWallView(weWall(:,1)==wNodeInd) = 1;
        weWallView(weWall(:,1)==eNodeInd) = -1;
    end
    
    % south and north wall
    [a, b, c]=meshgrid(wNodeInd:eNodeInd-1,[sNodeInd nNodeInd],tNodeInd:bNodeInd-1);
    snWall = [a(:) b(:) c(:)];
    if sNodeInd == nNodeInd % for cross section, sign determined by the inquiry location relative to the node
        snWallView = a(:)*0 + 2 * sign(mean(block(p,3:4)) - nodeY(sNodeInd) + 1e-10); % in case inquiry location exactly on the node, give a slight positive push so south wall dominates
    else
        snWallView = a(:)*0;
        snWallView(snWall(:,2)==sNodeInd) = 2;
        snWallView(snWall(:,2)==nNodeInd) = -2;
    end
    
    % top and bottom wall
    [a, b, c]=meshgrid(wNodeInd:eNodeInd-1,sNodeInd:nNodeInd-1,[tNodeInd bNodeInd]);
    tbWall = [a(:) b(:) c(:)];
    if tNodeInd == bNodeInd % for cross section, sign determined by the inquiry location relative to the node
        tbWallView = a(:)*0 + 3 * sign(mean(block(p,5:6)) - nodeZ(tNodeInd) - 1e-10); % in case inquiry location exactly on the node, give a slight negative push so top wall dominates
    else
        tbWallView = a(:)*0;
        tbWallView(tbWall(:,3)==tNodeInd) = -3;
        tbWallView(tbWall(:,3)==bNodeInd) = 3;
    end
  
    
    Wall{p} = [weWall; snWall; tbWall];
    WallView{p} = [weWallView; snWallView; tbWallView];
    
end


[a, b, c] = ndgrid(nodeZ,nodeX,nodeY);
node = [b(:) c(:) a(:) ]; % all nodes in a 3D mesh ordered using UBC convention

% how node index changes at vertex of a cell
% table value in relative increment
% vertex index   1 2   3     4       5         6           7
indIncreTable = [0 1 Nnz Nnz+1 Nnz*Nnx Nnz*Nnx+1 Nnz*Nnx+Nnz]';

vertexRule = [1 2 6 5;   % face orient type 1, normal direction x
              1 2 4 3;   % face orient type 2, normal direction y
              1 3 7 5];   % face orient type 3, normal direction z
          
figure;          
for p = 1:Nblk
    
    xyzind = Wall{p};
    orient = WallView{p};
          
    Nface = size(xyzind,1);
    nodeInd = zeros(Nface,4);
    
    % four vertices of face
    nodeInd(:,1) = Nnz*Nnx*(xyzind(:,2)-1) + Nnz*(xyzind(:,1)-1) + xyzind(:,3);
    nodeInd(:,2) = nodeInd(:,1) + indIncreTable(vertexRule(abs(orient),2));
    nodeInd(:,3) = nodeInd(:,1) + indIncreTable(vertexRule(abs(orient),3));
    nodeInd(:,4) = nodeInd(:,1) + indIncreTable(vertexRule(abs(orient),4));
    
    % find cell index using the sign of orient
    countCellY = xyzind(:,2) - 1 - (orient==-2);
    countCellX = xyzind(:,1) - 1 - (orient==-1);
    countCellZ = xyzind(:,3)     - (orient== 3);
    
    % safeguard beyond-boundary index
    countCellX(countCellX<0) = 0; countCellX(countCellX>Ncx-1) = Ncx-1;
    countCellY(countCellY<0) = 0; countCellY(countCellY>Ncy-1) = Ncy-1;
    countCellZ(countCellZ==0) = 1; countCellZ(countCellZ>Ncz) = Ncz;
    
    
    
    cellInd = Ncx * Ncz * countCellY + Ncz * countCellX + countCellZ;
    cdata = model(cellInd);
    
    hp(p) = patch('Faces',nodeInd,'Vertices',node,'FaceVertexCData',cdata,...
        'CDataMapping','scaled','FaceColor','flat');

end

% set(gca,'box','on','projection','perspective','DataAspectRatio',[1 1 1]);
set(gca,'box','on','DataAspectRatio',[1 1 1]);
grid on;
view([0.5 -1 1]);
hc = colorbar;

% make log-scale color bar
% if strcmpi(logflag,'log')
%     textexp = @(x)strcat('10^{',x,'}');
%     hc.TickLabels = textexp(hc.TickLabels);
% end


end