% Add round disks to model
% FUNCTION model_out = addDisk(nodeX,nodeY,nodeZ,model_in,dskLocInfo,dskValInfo)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter
%     model_in: input model vector, to which new blocks are added, if omitted, assign
%     0 everywhere in the mesh; can be a scalar
%     dskLocInfo: disk location info, Ndsk x 6 matrix, for each row 
%     [center_x center_y center_z radius thickness axis_direction];
%     axis_direction: 1 = x, 2 = y, 3 = z
%     dskValInfo: disk's model value, Ndsk x 1 vector, each element for one conductivity
% OUTPUT
%     model_out: model vector with objects embedded
% LAST MODIFIED 20191107 yangdikun@gmail.com
function model_out = addDisk(nodeX,nodeY,nodeZ,model_in,dskLocInfo,dskValInfo)

[Nx, Ny, Nz, x0, y0, z0, hx, hy, hz, centerX, centerY, centerZ] ...
    = getMeshPara(nodeX,nodeY,nodeZ);

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

Ndsk = size(dskLocInfo,1);
xyz = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,1:Nx*Ny*Nz);
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
    
for i = 1:Ndsk
    
    centerx = dskLocInfo(i,1);
    centery = dskLocInfo(i,2);
    centerz = dskLocInfo(i,3);
    radius = dskLocInfo(i,4);
    thick = dskLocInfo(i,5);
    direct = dskLocInfo(i,6);
    
    switch direct
        case 1 % axis pointing x
            ind = centerx-thick/2 <= x & x<=centerx+thick/2 & ...
                  sqrt((y-centery).^2+(z-centerz).^2)<=radius;
        case 2 % axis pointing y
            ind = centery-thick/2 <= y & y<=centery+thick/2 & ...
                  sqrt((x-centerx).^2+(z-centerz).^2)<=radius;
        case 3 % axis pointing z
            ind = centerz-thick/2 <= z & z<=centerz+thick/2 & ...
                  sqrt((x-centerx).^2+(y-centery).^2)<=radius;
    end

    model_out(ind) = dskValInfo(i);
    
end


end






