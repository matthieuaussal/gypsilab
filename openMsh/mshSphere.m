function mesh = mshSphere(N,rad)
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
%| Synopsis   : Build uniform mesh for a sphere                           |
%+========================================================================+

% Fibonnacci rules
or      = (1+sqrt(5))/2;
theta   = (mod((2*pi/or) .* (0:N-1),2*pi))';
phi     = asin( -1 + 2/(N-1) * (0:N-1))';

% Carthesian coordinates
[x,y,z] = sph2cart(theta,phi,rad);
X       = [0 0 0 ; x y z];

% Delaunay triangulation
DT        = delaunayTriangulation(X);
[elt,vtx] = freeBoundary(DT);

% Mesh
mesh = msh(vtx,elt);
end
