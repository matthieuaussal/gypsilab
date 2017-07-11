function [X,W,elt2qud] = mshQuadrature(mesh)
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
%|              François Alouges - CMAP, Ecole Polytechnique              |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Apply reference quadrature to any mesh                    |
%+========================================================================+
    
% Reference quadrature
[x,w] = mshReference(mesh);

% Initialization
X       = zeros(size(x,1)*size(mesh.elt,1),size(mesh.vtx,2));
elt2qud = zeros(size(mesh.elt,1),size(x,1));

% For each quadrature point
for j = 1:size(x,1)
    % Indice
    idx          = (j:size(x,1):size(x,1)*size(mesh.elt,1))';
    elt2qud(:,j) = idx;
    
    % Particles mesh
    if (size(mesh.elt,2) == 1)
        error('mshQuadrature.m : unavailable case')
            
    % Edge mesh
    elseif (size(mesh.elt,2) == 2)
        X(idx,:) = (1-x(j,1)) * mesh.vtx(mesh.elt(:,1),:) ...
            + x(j,1) * mesh.vtx(mesh.elt(:,2),:) ;

    % Triangular mesh
    elseif (size(mesh.elt,2) == 3)
        X(idx,:) = (1-x(j,1)-x(j,2)) * mesh.vtx(mesh.elt(:,1),:) ...
            + x(j,1) * mesh.vtx(mesh.elt(:,2),:) ...
            + x(j,2) * mesh.vtx(mesh.elt(:,3),:);
        
    % Tetrahedron mesh
    elseif (size(mesh.elt,2) == 4)
        X(idx,:) = (1-x(j,1)-x(j,2)-x(j,3)) * mesh.vtx(mesh.elt(:,1),:) ...
            + x(j,1) * mesh.vtx(mesh.elt(:,2),:) ...
            + x(j,2) * mesh.vtx(mesh.elt(:,3),:) ...
            + x(j,3) * mesh.vtx(mesh.elt(:,4),:);
        
    % Unknown type
    else
        error('mshQuadrature.m : unavailable case')
    end
end

% Quadrature weight
W = kron(mesh.ndv,w);
end