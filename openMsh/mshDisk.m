function mesh = mshDisk(N,rad)
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
%| Synopsis   : Build uniform mesh for a disk                             |
%+========================================================================+

% Radial discretisation
dr = sqrt(pi*rad^2/N);
dr = rad/ceil(rad/dr);
r  = dr:dr:rad;

% Angular uniform discretization
rho = cell(length(r),1); theta = rho;
for ir = 1:length(r)
    dtheta = dr/r(ir);
    dtheta = 2*pi/ceil(2*pi/dtheta);    
    theta{ir} = (0:dtheta:2*pi-dtheta)';
    rho{ir}   = r(ir)*ones(length(theta{ir}),1);
end

% Carthesian coordinates
[x,y] = pol2cart(cell2mat(theta),cell2mat(rho));
X     = [0 0 ; x y];

% Unicity test
tmp = unique(X,'rows','stable');
if max(abs(X-tmp)) > 1e-12
    error('mshDisk : non unicity of vertices')
end
   
% Delaunay triangulation
DT = delaunayTriangulation(x,y);

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(vtx,elt);
end
