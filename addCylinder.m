function model_out = addCylinder(nodeX, nodeY, nodeZ, model_in, cldLocInfo, cldValInfo)
% Add round cylinders to model
% FUNCTION model_out = addCylinder(nodeX,nodeY,nodeZ,model_in,cldLocInfo,cldValInfo)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter
%     model_in: input model vector, to which new blocks are added, if omitted, %     assign 0 everywhere in the mesh; can be a scalar
%     cldLocInfo: cylinder location info, Ncld x 7 matrix, for each row 
%     [endpoint_x1 end_point_y1 endpoint_z1 endpoint_x2 endpoint_y2 endpoint_z2 radius];
%     endpoint(_x1, y1, z1, x2, y2, z2): end point coordinates of cylinderical axis
%     dcldValInfo: cylinder's model value, Ncld x 1 vector, each element for one conductivity
% OUTPUT
%     model_out: model vector with objects embedded
% LAST MODIFIED 20201201 ycli0536@gmail.com

Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

Ncld = size(cldLocInfo,1);
temp = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,[]);
x = temp(:,1);
y = temp(:,2);
z = temp(:,3);

syms a b c

% Rotataion matrix (x-y-z)
R = [ cos(b)*cos(c)                       cos(b)*sin(c)                      -sin(b);
     -cos(a)*sin(c)+sin(a)*sin(b)*cos(c)  cos(a)*cos(c)+sin(a)*sin(b)*sin(c)  sin(a)*cos(b);
      sin(a)*sin(c)+cos(a)*sin(b)*cos(c) -sin(a)*cos(c)+cos(a)*sin(b)*sin(c)  cos(a)*cos(b)];

for i = 1:Ncld

    x1 = cldLocInfo(i,1);
    y1 = cldLocInfo(i,2);
    z1 = cldLocInfo(i,3);
    x2 = cldLocInfo(i,4);
    y2 = cldLocInfo(i,5);
    z2 = cldLocInfo(i,6);
    radius = cldLocInfo(i,7);
    
    cylinder_axis = [x2 - x1; y2 - y1; z2 - z1];
    direct = cylinder_axis / norm(cylinder_axis);
    f = R * [0; 0; 1] - direct;
    func = matlabFunction(f, 'Vars',{[a, b, c]});
    solution = fsolve(func, [0, 0, 0]);
    
    alpha = solution(1);
    beta = solution(2);
    theta = solution(3);
    Len = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    % displacement_vector: (mid point of the cylinder axis - Origin)
    displacement_vector = [(x1+x2)/2 (y1+y2)/2 (z1+z2)/2];
    
    x0 = displacement_vector(1);
    y0 = displacement_vector(2);
    z0 = displacement_vector(3);
    
    xt = - x0*cos(beta)*cos(theta)  + x*cos(beta)*cos(theta) ...
         + y0*cos(alpha)*sin(theta) - y*cos(alpha)*sin(theta) ...
         - z0*sin(alpha)*sin(theta) + z*sin(alpha)*sin(theta) ...
         - z0*cos(alpha)*sin(beta)*cos(theta) + z*cos(alpha)*sin(beta)*cos(theta) ...
         - y0*sin(alpha)*sin(beta)*cos(theta) + y*sin(alpha)*sin(beta)*cos(theta);
    yt = - y0*cos(alpha)*cos(theta) + y*cos(alpha)*cos(theta)...
         - x0*cos(beta)*sin(theta)  + x*cos(beta)*sin(theta) ...
         + z0*sin(alpha)*cos(theta) - z*sin(alpha)*cos(theta) ...
         - z0*cos(alpha)*sin(beta)*sin(theta) + z*cos(alpha)*sin(beta)*sin(theta) ...
         - y0*sin(alpha)*sin(beta)*sin(theta) + y*sin(alpha)*sin(beta)*sin(theta);
    zt = - x*sin(beta) + y*cos(beta)*sin(alpha) + z*cos(alpha)*cos(beta) ...
         + x0*sin(beta) - z0*cos(alpha)*cos(beta) - y0*cos(beta)*sin(alpha);
    
    % Assumed standard cylinder expression: (x/a)^2 + (y/b)^2 = 1, where z
    % from [-Z, Z]
    dist = sqrt( (xt.^2)./radius^2 + (yt.^2)./radius^2);
    ind = (dist <= 1) & (zt >= -Len/2) & (zt <= Len/2);
    
    model_out(ind) = cldValInfo(i);

end

end