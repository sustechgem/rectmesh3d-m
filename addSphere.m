% Add a sphere to model
% FUNCTION model_out = addSphere(nodeX,nodeY,nodeZ,model_in,sphereLocInfo,sphereValInfo)
% INPUT
%     nodeX,nodeY,nodeZ,model_in: input mesh and model
%     sphereLocInfo: sphere location info, Nsphere x 5 matrix, for each row [x,y,z,radius_in,radius_out]
%       center: center [x,y,z] of sphere
%       radius_in, radius_out: inner and outer radii of sphere (upper and
%       lower bounds of transition zone)
%     sphereValInfo: sphere's value, Nblk x 2 vector, for each row [value_in,value_out]
%       val_in, val_out: model values for the sphere (in) and the outmost 
%       edge of transition zone (out); cells outside of radius_out will not be changed
% OUTPUT
%     model_out: model vector with objects embedded
% NOTE
%     Smartly use the inner and outer radius to make a smooth-edge sohere
%     or a spherical shell
% LAST MODIFIED 20191107 yangdikun@gmail.com
function model_out = addSphere(nodeX,nodeY,nodeZ,model_in,sphereLocInfo,sphereValInfo)

Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

Nsphere = size(sphereLocInfo,1);

for p = 1:Nsphere
    
    center = sphereLocInfo(p,1:3);
    radius_in = sphereLocInfo(p,4);
    radius_out = sphereLocInfo(p,5);
    val_in = sphereValInfo(p,1);
    val_out = sphereValInfo(p,2);
    
    temp = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,[]);
    x = temp(:,1);
    y = temp(:,2);
    z = temp(:,3);
    
    dist = sqrt( (x-center(1)).^2 + (y-center(2)).^2 + (z-center(3)).^2);
    
    model_out(dist<=radius_in) = val_in;
    
    f = @(x)(exp(-1./x));
    g = @(x,a,b)(1- f((x-a)./(b-a)) ./ ( f((x-a)./(b-a)) + f(1-(x-a)./(b-a))  ) );
    
    ind = dist>radius_in & dist<radius_out; % transition part
    model_out(ind) = g(dist(ind),radius_in,radius_out) * (val_in-val_out) + val_out;
end
    
end