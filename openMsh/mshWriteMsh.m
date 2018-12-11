function mshWriteMsh(filename,mesh)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshWriteMsh.m                                 |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Write mesh and data to msh format             |
%|  `---'  |                (particle, edge, triangular and tetrahedral)  |
%+========================================================================+

% Security
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) > 4)
    error('mshWriteMsh.m : unavailable case')
end
   
% Open file
fid = fopen(filename,'w');

% Header
fprintf(fid,'%s\n','$MeshFormat');
fprintf(fid,'%s\n','2.2 0 8');
fprintf(fid,'%s\n','$EndMeshFormat');

% Nodes
fprintf(fid,'%s\n','$Nodes');
fprintf(fid,'%d\n',size(mesh.vtx,1));
for i = 1:size(mesh.vtx,1)
    fprintf(fid,'%d %f %f %f\n',i,mesh.vtx(i,:));
end
fprintf(fid,'%s\n','$EndNodes');

% Elements (up to tetra)
fprintf(fid,'%s\n','$Elements');
fprintf(fid,'%d\n',size(mesh.elt,1));
for i = 1:size(mesh.elt,1)
    if (size(mesh.elt,2) == 1)  % particle
        fprintf(fid,'%d %d %d %d %d  %d\n',i,15,2,0,1,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 2)  % segment
        fprintf(fid,'%d %d %d %d %d  %d %d\n',i,1,2,0,1,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 3)  % triangle
        fprintf(fid,'%d %d %d %d %d  %d %d %d\n',i,2,2,0,1,mesh.elt(i,:));
    elseif (size(mesh.elt,2) == 4)  % tetra
        fprintf(fid,'%d %d %d %d %d  %d %d %d %d\n',i,4,2,0,1,mesh.elt(i,:));
    end
end
fprintf(fid,'%s\n','$EndElements');

% Close file
fclose(fid);
disp([filename,' created.']);
end
