% Make a 3D rectilinear mesh; 
% Mesh is of core plus padding structure and adjustable in all 6 directions
% FUNCTION [nodeX,nodeY,nodeZ,topomnz] = makeRectMesh3D(origin,Wsize0,Wrate,Wpadding,Esize0,Erate,Epadding,...
%    Ssize0,Srate,Spadding,Nsize0,Nrate,Npadding,Tsize0,Trate,Tpadding,Bsize0,Brate,Bpadding,topoxyz)
% INPUT
%     origin: [x y z] of the "center" (origin) of the mesh
%     ?size0, ?rate, ?padding: smallest cell size, expansion rate and padding
%     distance at west/east/south/north/top/bottom directions
%     topoxyz: 3-column xyz matrix for the topography (active or inactive
%     cells); if empty, use the z coordinate of the origin for a flat topo
% OUTPUT
%     nodeX,nodeY,nodeZ: vectors of nodes in x, y, z that define a 3D mesh
%     topomnz: a vector for topo (active cell = 1; inactive cell = 0)
% LAST MODIFIED 20191122 yangdikun@gmail.com
function [nodeX,nodeY,nodeZ,topomnz] = makeRectMesh3D(origin,Wsize0,Wrate,Wpadding,Esize0,Erate,Epadding,...
    Ssize0,Srate,Spadding,Nsize0,Nrate,Npadding,Tsize0,Trate,Tpadding,Bsize0,Brate,Bpadding,topoxyz)


% west
W = reshape(Wsize0,[],1);
while sum(W) < Wpadding + sum(Wsize0)
    W = [W; W(end)*Wrate];
end
% east
E = reshape(Esize0,[],1);
while sum(E) < Epadding + sum(Esize0)
    E = [E; E(end)*Erate];
end
hx = [flipud(W); E];
x0 = origin(1) - sum(W);
nodeX = size2node(x0,hx,'x');

% south
S = reshape(Ssize0,[],1);
while sum(S) < Spadding + sum(Ssize0)
    S = [S; S(end)*Srate];
end
% north
N = reshape(Nsize0,[],1);
while sum(N) < Npadding + sum(Nsize0)
    N = [N; N(end)*Nrate];
end
hy = [flipud(S); N];
y0 = origin(2) - sum(S);
nodeY = size2node(y0,hy,'y');

% top
T = reshape(Tsize0,[],1);
while sum(T) < Tpadding + sum(Tsize0)
    T = [T; T(end)*Trate];
end
% east
B = reshape(Bsize0,[],1);
while sum(B) < Bpadding + sum(Bsize0)
    B = [B; B(end)*Brate];
end
hz = [flipud(T); B];
z0 = origin(3) + sum(T);
nodeZ = size2node(z0,hz,'z');

% make topo
if isempty(topoxyz) % constant topo = elev
    topomnz = repmat([ T*0; B*0+1 ],(length(nodeX)-1)*(length(nodeY)-1),1);
else % 3-column topo
    ind = topoxyz(:,1)>=nodeX(1) & topoxyz(:,1)<=nodeX(end) ...
        & topoxyz(:,2)>=nodeY(1) & topoxyz(:,2)<=nodeY(end);
    topomnz = PointTopo2MeshTopo(topoxyz(ind,:),nodeX,nodeY,nodeZ);
end


end