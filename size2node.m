% Calculate cell node in one direction of a 3D mesh
% FUNCTION node = size2node(x0,hx,flag)
% INPUT
%     x0: coordinate of the first node
%     hx: a vector of cell sizes
%     flag: 'x', 'y', or 'z' to indicate the direction
% OUTPUT
%     node: a vector of node in a particular direction +x/+y/-z
% LAST MODIFIED 20191122 yangdikun@gmail.com
function node = size2node(x0,hx,flag)
    hx = reshape(hx,length(hx),1);
    if strcmpi(flag,'z')
        node = [x0; x0-cumsum(hx)];
    else
        node = [x0; x0+cumsum(hx)];
    end
end
        

