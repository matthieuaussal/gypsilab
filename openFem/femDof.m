function [X,elt2dof] = femDof(fe,mesh)
%+========================================================================+
%|                                                                        |
%|                  OPENFEM, FINITE AND BOUNDARY ELEMENT                  |
%|              openFem is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openFem,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%|              François Alouges - CMAP, Ecole Polytechnique              |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Define degrees of freedom                                 |
%+========================================================================+

% Lagrange order 0, constant by element    
if strcmp(fe.typ,'P0')
    X       = mesh.ctr;
    elt2dof = (1:size(mesh.elt,1))';

% Lagrange order 1, piecewise linear by element    
elseif strcmp(fe.typ,'P1')
    X       = mesh.vtx;
    elt2dof = mesh.elt;
    
% Lagrange order 2, piecewise quadratic by element    
elseif strcmp(fe.typ,'P2')
    [edg2vtx,elt2edg] = mesh.edg;
    X                 = [mesh.vtx ; 0.5.*(mesh.vtx(edg2vtx(:,1),:)+mesh.vtx(edg2vtx(:,2),:))];
    elt2dof           = [mesh.elt , elt2edg + size(mesh.vtx,1)];

% Others    
else
    error('femDof.m : unavailable case')
end
end
