% Add Ellipsoid to model
% FUNCTION model_out = addEllipsoid(nodeX,nodeY,nodeZ,model_in,ellipLocInfo,ellipValInfo)
% INPUT
%     nodeX,nodeY,nodeZ: mesh parameter
%     model_in: input model, to which new blocks are added, if omitted, assign
%     ellipLocInfo: ellipsoid location info, Nellip x 9 matrix,[x,y,z,alpha,beta,theta,a,b,c]
%       [x, y, z]: the center of ellipsoid,
%       [alpha, beta, theta]: the angle of ellipsoid, angle range is [-pi/2,pi/2]
%       [a, b, c]: the lenth of three axis
%     ellipValInfo: ellipsoid's value, Nblk x 1 vector, each element for one block
% OUTPUT
%     model_out: model vector with objects embedded
% LAST MODIFIED 20210910 Huying.em@gmail.com
function model_out = addEllipsoid(nodeX,nodeY,nodeZ,model_in,ellipLocInfo,ellipValInfo)
Nx = length(nodeX) - 1;
Ny = length(nodeY) - 1;
Nz = length(nodeZ) - 1;

if isempty(model_in)
    model_in = zeros(Nx*Ny*Nz,1);
elseif isscalar(model_in)
    model_in = zeros(Nx*Ny*Nz,1) + model_in;
end
model_out = model_in;

Nellip = size(ellipLocInfo,1);

for p = 1:Nellip
    
    center = ellipLocInfo(p,1:3);
    angle = ellipLocInfo(p,4:6);
    axis = ellipLocInfo(p,7:9);
    val_out = ellipValInfo(p);
    
    temp = CellIndex2PointXYZ(nodeX,nodeY,nodeZ,[]);
    x = temp(:,1);
    y = temp(:,2);
    z = temp(:,3);
    
    x0 = center(1);
    y0 = center(2);
    z0 = center(3);
    
    alpha = angle(1);
    beta = angle(2);
    theta = angle(3);
    
    radius_a = axis(1);
    radius_b = axis(2);
    radius_c = axis(3);
    
    x_p = (z*cos(alpha)*sin(beta)*cos(theta) - z0*cos(alpha)*sin(beta)*cos(theta) + y*sin(alpha)*sin(beta)*cos(theta)...
        - y0*sin(alpha)*sin(beta)*cos(theta) + x*cos(alpha)^2*cos(beta)*cos(theta) - x0*cos(alpha)^2*cos(beta)*cos(theta) ...
        + x*cos(beta)*sin(alpha)^2*cos(theta) - x0*cos(beta)*sin(alpha)^2*cos(theta) - y*cos(alpha)*cos(beta)^2*sin(theta)...
        + y0*cos(alpha)*cos(beta)^2*sin(theta) - y*cos(alpha)*sin(beta)^2*sin(theta) + y0*cos(alpha)*sin(beta)^2*sin(theta) ...
        + z*cos(beta)^2*sin(alpha)*sin(theta) - z0*cos(beta)^2*sin(alpha)*sin(theta) + z*sin(alpha)*sin(beta)^2*sin(theta) - ....
        z0*sin(alpha)*sin(beta)^2*sin(theta))/((cos(theta)^2 + sin(theta)^2)*(cos(alpha)^2*cos(beta)^2 + cos(alpha)^2*sin(beta)^2 ...
        + cos(beta)^2*sin(alpha)^2 + sin(alpha)^2*sin(beta)^2));
    y_p = (z*cos(alpha)*sin(beta)*sin(theta) - z0*cos(alpha)*sin(beta)*sin(theta) + y*sin(alpha)*sin(beta)*sin(theta)...
        - y0*sin(alpha)*sin(beta)*sin(theta) + y*cos(alpha)*cos(beta)^2*cos(theta) - y0*cos(alpha)*cos(beta)^2*cos(theta)...
        + x*cos(alpha)^2*cos(beta)*sin(theta) - x0*cos(alpha)^2*cos(beta)*sin(theta) + y*cos(alpha)*sin(beta)^2*cos(theta)...
        - y0*cos(alpha)*sin(beta)^2*cos(theta) - z*cos(beta)^2*sin(alpha)*cos(theta) + z0*cos(beta)^2*sin(alpha)*cos(theta)...
        + x*cos(beta)*sin(alpha)^2*sin(theta) - x0*cos(beta)*sin(alpha)^2*sin(theta) - z*sin(alpha)*sin(beta)^2*cos(theta)...
        + z0*sin(alpha)*sin(beta)^2*cos(theta))/((cos(theta)^2 + sin(theta)^2)*(cos(alpha)^2*cos(beta)^2 + cos(alpha)^2*sin(beta)^2 ...
        + cos(beta)^2*sin(alpha)^2 + sin(alpha)^2*sin(beta)^2));
    z_p = (z*cos(alpha)*cos(beta) - z0*cos(alpha)*cos(beta) + y*cos(beta)*sin(alpha) - y0*cos(beta)*sin(alpha)...
        - x*cos(alpha)^2*sin(beta) + x0*cos(alpha)^2*sin(beta) - x*sin(alpha)^2*sin(beta) + x0*sin(alpha)^2*sin(beta))...
        /(cos(alpha)^2*cos(beta)^2 + cos(alpha)^2*sin(beta)^2 + cos(beta)^2*sin(alpha)^2 + sin(alpha)^2*sin(beta)^2);
    
    dist = sqrt( (x_p.^2)./radius_a^2 + (y_p.^2)./radius_b^2 + (z_p.^2)./radius_c^2 );
    ind = dist<=1;
    model_out(ind) = val_out;
end

end