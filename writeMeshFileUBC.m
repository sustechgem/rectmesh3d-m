% write UBC-GIF mesh file to disk
% FUNCTION writeMeshFileUBC(meshName,nodeX,nodeY,nodeZ)
% INPUT
%     meshFileName: mesh file name in a string (text file in UC-GIF format)
%     nodeX, nodeY, nodeZ: mesh description in nodes
% LAST MODIFIED 20210908 yangdikun@gmail.com
function writeMeshFileUBC(meshName,nodeX,nodeY,nodeZ)

[Nx, Ny, Nz, x0, y0, z0, hx, hy, hz] = getMeshPara(nodeX,nodeY,nodeZ);

fid = fopen(meshName,'wt');

fprintf(fid,'%d %d %d\n',Nx,Ny,Nz);
fprintf(fid,'%f %f %f\n',x0,y0,z0);
fprintf(fid,'%f ',hx);
fprintf(fid,'\n');
fprintf(fid,'%f ',hy);
fprintf(fid,'\n');
fprintf(fid,'%f ',hz);
fprintf(fid,'\n');
fclose(fid);

end




