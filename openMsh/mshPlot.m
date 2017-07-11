function mshPlot(data)
%+========================================================================+
%|                                                                        |
%|           OPENMSH, MESH MANAGEMENT AND NUMERICAL QUADRATURE            |
%|              openMsh is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openMsh,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Plot mesh with or without data                            |
%+========================================================================+

% Mesh 
mesh = data{1};

% Particles mesh
if (size(mesh.elt,2) == 1)
    x = mesh.vtx(mesh.elt,1);
    y = mesh.vtx(mesh.elt,2);
    z = mesh.vtx(mesh.elt,3);
    if (length(data)==1)
        plot3(x,y,z,'.k');
    else
        plot3(x,y,z,data{2});
    end

% Edge mesh
elseif (size(mesh.elt,2) == 2)
    x = [mesh.vtx(mesh.elt(:,1),1) mesh.vtx(mesh.elt(:,2),1)]';
    y = [mesh.vtx(mesh.elt(:,1),2) mesh.vtx(mesh.elt(:,2),2)]';
    z = [mesh.vtx(mesh.elt(:,1),3) mesh.vtx(mesh.elt(:,2),3)]';
    if (length(data)==1)
        plot3(x,y,z,'.-k');
    else
        plot3(x,y,z,data{2})
    end
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    if (length(data) == 1)
        V = zeros(size(mesh.elt,1),1);
    elseif (numel(data{2}) == 1)
        V = data{2} * ones(size(mesh.elt,1),1);
    else
        V = data{2};
    end
    h = trisurf(mesh.elt,mesh.vtx(:,1),mesh.vtx(:,2),mesh.vtx(:,3),V);
    if (length(data) == 1)
        set(h,'EdgeColor','k');
        set(h,'FaceColor','w');
        hold on
        plot3(mesh.vtx(:,1),mesh.vtx(:,2),mesh.vtx(:,3),'.k')
        hold off
    else
        set(h,'EdgeColor','none');
    end   
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    fce2vtx = mesh.fce;
    data{1} = msh(mesh.vtx,fce2vtx);
    mshPlot(data);
    
% Unknown type    
else
    error('mshPlot.m : unavailable case')
end
end
