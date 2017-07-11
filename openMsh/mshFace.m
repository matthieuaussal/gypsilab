function [fce2vtx,elt2fce] = mshFace(mesh)
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
%| Synopsis   : Compute face table                                        |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshFace.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    error('mshFace.m : unavailable case')
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    fce2vtx = mesh.elt;
    elt2fce = (1:size(mesh.elt,1))';
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    % All faces
    fce2vtx = [ mesh.elt(:,[1,2,3]) ; mesh.elt(:,[2,3,4]) ; ...
        mesh.elt(:,[3,4,1]) ; mesh.elt(:,[4,1,2]) ];
    
    % Faces unicity
    tmp           = sort(fce2vtx,2);
    [~,I,elt2fce] = unique(tmp,'rows');
    
    % Final faces
    fce2vtx = fce2vtx(I,:);
    elt2fce = reshape(elt2fce,size(mesh.elt,1),size(mesh.elt,2));
    
% Unknown type
else
    error('mshEdge.m : unavailable case')
end

end
