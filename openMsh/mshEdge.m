function [edg2vtx,elt2edg] = mshEdge(mesh)
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
%| Synopsis   : Compute edges tables                                      |
%+========================================================================+

% Particles mesh
if (size(mesh.elt,2) == 1)
    error('mshEdge.m : unavailable case')
    
% Edge mesh
elseif (size(mesh.elt,2) == 2)
    edg2vtx = mesh.elt;
    elt2edg = (1:size(mesh.elt,1))';
    
% Triangular mesh
elseif (size(mesh.elt,2) == 3)
    % All edges
    edg2vtx = [ mesh.elt(:,[1 2]) ; mesh.elt(:,[2 3]) ; mesh.elt(:,[3 1]) ];
    
    % Edges unicity
    tmp           = sort(edg2vtx,2);
    [~,I,elt2edg] = unique(tmp,'rows');
    
    % Final Edges
    edg2vtx = edg2vtx(I,:);
    elt2edg = reshape(elt2edg,size(mesh.elt,1),size(mesh.elt,2));
    
% Tetrahedron mesh
elseif (size(mesh.elt,2) == 4)
    % All edges
    edg2vtx = [ mesh.elt(:,[1 2]) ; mesh.elt(:,[2 3]) ; mesh.elt(:,[3 1]) ; 
        mesh.elt(:,[1 4]) ; mesh.elt(:,[2 4]) ; mesh.elt(:,[3 4]) ]; 
        
    % Edges sunicity
    tmp           = sort(edg2vtx,2);
    [~,I,elt2edg] = unique(tmp,'rows');
    
    % Final Edges
    edg2vtx = edg2vtx(I,:);
    elt2edg = reshape(elt2edg,size(mesh.elt,1),6);
    
% Unknown type
else
    error('mshEdge.m : unavailable case')
end

end
