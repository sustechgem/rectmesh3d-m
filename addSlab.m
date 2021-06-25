% Add a Slab to model
% FUNCTION model_out = addSlab(nodeX,nodeY,nodeZ,model_in,sfLocInfo,th, val)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter nodeX ¡Ê[Xmin, Xmax]£¬nodeY ¡Ê[Ymin, Ymax]
%     model_in: input model, to which new blocks are added, if omitted, assign
%     0 everywhere in the mesh; can be a scalar
%     sfLocInfo: surface location info, a matrix containing the coordinates
%     [x y z] of a set of scattered points about a surface representing a
%     plane at half the thickness of the slab
%     th: the thickness of the slab
%     val: model value to be assgined to the target cells
% OUTPUT
%     model_out: model vector with objects embedded
% NOTE
%     The surface represented by sfLocInfo cannot be parallel to the XOZ and
%     YOZ plane
% LAST MODIFIED 20210528 yinchu.li@hotmail.com
function model_out = addSlab(nodeX, nodeY, nodeZ, model_in, sfLocInfo, th, val)
    Nx = length(nodeX) - 1;
    Ny = length(nodeY) - 1;
    Nz = length(nodeZ) - 1;

    if isempty(model_in)
        model_in = zeros(Nx*Ny*Nz,1);
    elseif isscalar(model_in)
        model_in = zeros(Nx*Ny*Nz,1) + model_in;
    end
    model_out = model_in;

    sfLocInfo_upper = sfLocInfo;
    sfLocInfo_bottom = sfLocInfo;
    sfLocInfo_upper(:, 3) = sfLocInfo(:, 3) + th / 2;
    sfLocInfo_bottom(:, 3) = sfLocInfo(:, 3) - th / 2;
    model_upperSurf = addSurface(nodeX,nodeY,nodeZ,model_in,sfLocInfo_upper,val);
    model_bottomSurf = addSurface(nodeX,nodeY,nodeZ,model_in,sfLocInfo_bottom,val);
    ind = (model_upperSurf - model_bottomSurf ~= 0);
    model_out(ind) = val;
end
