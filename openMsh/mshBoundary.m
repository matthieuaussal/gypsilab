function bnd2vtx = mshBoundary(mesh)
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
%| Synopsis   : Extract element from boundary                             |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshBoundary.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    edg2vtx = mesh.edg;
    mlt     = accumarray(edg2vtx(:),1,[size(mesh.vtx,1),1]);
    bnd2vtx = find(mlt==1);
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    dt      = triangulation(mesh.elt,mesh.vtx);
    bnd2vtx = freeBoundary(dt);
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    dt      = triangulation(mesh.elt,mesh.vtx);
    bnd2vtx = freeBoundary(dt);
    
else
    error('mshBoundary.m : unavailable case')
end
end
