function mesh = mshSquare(N,L)
%+========================================================================+
%|                                                                        |
%|           OPENMSH, MESH MANAGEMENT AND NUMERICAL QUADRATURE            |
%|              openMsh is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
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
%|              François Alouges - CMAP, Ecole polytechnique              |
%| Creation   : 06.07.17                                                  |
%| Last modif : 06.07.17                                                  |
%| Synopsis   : Build uniform mesh for a unit square                      |
%+========================================================================+

% Optimal number of point for each dimension
n = L/min(L);
n = round( n * (N/prod(n))^(1/2) );

% Delaunay mesh
x = -0.5*L(1) : L(1)/n(1) : 0.5*L(1);
y = -0.5*L(2) : L(2)/n(2) : 0.5*L(2);

% Triangulation
[x,y]  = meshgrid(x,y);
z      = zeros(size(x));
DT     = delaunayTriangulation([x(:) y(:)]);
Points = [DT.Points z(:)];

% Build mesh
mesh = msh(Points,DT.ConnectivityList);
end