% Add a polyhedron to model
% FUNCTION model_out = addPolyhedron(nodeX,nodeY,nodeZ,model_in,points,val)
% INPUT
%     nodeX,nodeY,nodeZ,model_in: input mesh and model
%     points: a 3-column matrix to specify a group of X-Y-Z points that 
%     form a convex hull enclosing target cells' centers
%     val: model value to be assgined to the target cells
% OUTPUT
%     model_out: model vector with objects embedded
% NOTE
%     Only one object at one function call
% LAST MODIFIED 20191107 yangdikun@gmail.com
function model_out = addPolyhedron(nodeX,nodeY,nodeZ,model_in,points,val)

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

tri = delaunayn(points);
tn = tsearchn(points,tri,[x y z]);
isInside = ~isnan(tn);

model_out(isInside) = val;

end