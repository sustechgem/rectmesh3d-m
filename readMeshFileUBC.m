% read UBC-GIF 3D mesh file from disk
% FUNCTION [nodeX, nodeY, nodeZ] = readMeshFileUBC(meshFileName)
% INPUT
%     meshFileName: mesh file name in a string (text file in UC-GIF format)
% OUTPUT
%     nodeX, nodeY, nodeZ: mesh description in nodes
% LAST MODIFIED 20210908 yangdikun@gmail.com
function [nodeX, nodeY, nodeZ] = readMeshFileUBC(meshFileName)

fid = fopen(meshFileName,'r');

temp = fgetl(fid);
v = sscanf(temp,'%f',3);
Nx = v(1);     %East 
Ny = v(2);     %North
Nz = v(3);     %Vertical

temp = fgetl(fid);
v = sscanf(temp,'%f',3);
x0 = v(1);
y0 = v(2);
z0 = v(3);

hx=[]; hy=[]; hz=[];

for j=1:3
    % Get line:
    t = fgetl(fid);
    % Mark end of line:
    t(end+1) = '!';
    go = 1;
    while go
        % Get number:
        [num,count,errmsg,nextindex] = sscanf(t,'%f',1);
        if ~isempty(errmsg) || count~=1
            fclose(fid);
            return
        end
        t = t(nextindex:end);
        t = strtrim(t);
        % Is it a single number:
        if t(1)=='*'
            % The number was the times:
            times = num;
            % Get number:
            t = t(2:end);
            t = strtrim(t);
            [num,count,errmsg,nextindex] = sscanf(t,'%f',1);
            if ~isempty(errmsg) || count~=1
                fclose(fid);
                return
            end
            t = t(nextindex:end);
            t = strtrim(t);
        else
            times = 1;
        end
        switch j
            case(1)
                hx = [hx ; repmat(num,times,1)];
            case(2)
                hy = [hy ; repmat(num,times,1)];
            case(3)
                hz = [hz ; repmat(num,times,1)];
        end
        % Check for end of line:
        if t(1)=='!'
            go = 0;
            break
        end
    end
end
fclose(fid);

nodeX = size2node(x0,hx,'x');
nodeY = size2node(y0,hy,'y');
nodeZ = size2node(z0,hz,'z');


end