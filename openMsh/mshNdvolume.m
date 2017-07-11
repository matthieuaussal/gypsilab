function S = mshNdvolume(mesh)
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
%| Synopsis   : Compute element volume of dimension n                     |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshSurface.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    
    % Surface
    S = sqrt(sum(E1.^2,2));
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
    E3 = cross(E1,E2) ;
    
    % Surface
    S = 0.5*sqrt(sum(E3.^2,2));
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    % Basis vector
    E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
    E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
    E3 = mesh.vtx(mesh.elt(:,4),:) - mesh.vtx(mesh.elt(:,1),:);

    % Surface
    S = 1/6 .* abs( ...
          E1(:,1) .* (E2(:,2).*E3(:,3) - E2(:,3).*E3(:,2)) ...
        - E2(:,1) .* (E1(:,2).*E3(:,3) - E1(:,3).*E3(:,2)) ...
        + E3(:,1) .* (E1(:,2).*E2(:,3) - E1(:,3).*E2(:,2)) );
    
% Unknown type
else
    error('mshSurface.m : unavailable case')
end

end
