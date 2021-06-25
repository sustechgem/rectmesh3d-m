% Add subsurface to model
% FUNCTION model_out = addSurface(nodeX,nodeY,nodeZ,model_in,sfLocInfo,val_out)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter nodeX ¡Ê[Xmin, Xmax]£¬nodeY ¡Ê[Ymin, Ymax]
%     model_in: input model, to which new blocks are added, if omitted, assign
%     0 everywhere in the mesh; can be a scalar
%     sfLocInfo: surface location info, a matrix containing the coordinates
%     [x y z] of a set of scattered points about a surface
%     val_out: model value under the surface (/layer)
% OUTPUT
%     model_out: model vector with objects embedded
% LAST MODIFIED 20210316 yinchu.li@hotmail.com
function model_out = addSurface(nodeX, nodeY, nodeZ, model_in, sfLocInfo, val_out)
Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

temp = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,[]);
x = temp(:,1);
y = temp(:,2);
z = temp(:,3);
xcenter = unique(x);
ycenter = unique(y);
zcenter = unique(z);

[xq,yq] = meshgrid(xcenter, ycenter);
zq = griddata(sfLocInfo(:, 1),sfLocInfo(:, 2),sfLocInfo(:,3),xq,yq);
zq_t = zq';

zLoc = repelem(zq_t(:), length(zcenter));
ind = z - zLoc <= 0;
model_out(ind) = val_out;

end