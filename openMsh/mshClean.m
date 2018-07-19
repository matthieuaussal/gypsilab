function mesh = mshClean(mesh,dst)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshClean.m                                    |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 21.06.2018                                    |
%| ( === ) |   SYNOPSIS   : Clean mesh                                    |
%|  `---'  |                                                              |
%+========================================================================+

% Unify duplicate vertex
if isempty(dst)
    [~,I,J] = unique(single(mesh.vtx),'rows','stable');
else
    [~,I,J] = unique(round(mesh.vtx/dst),'rows','stable');
end
mesh.vtx = mesh.vtx(I,:);
if (size(mesh.elt,1) == 1)
    J = J';
end
mesh.elt = J(mesh.elt);

% Extract vertex table from element
Ivtx              = zeros(size(mesh.vtx,1),1);
Ivtx(mesh.elt(:)) = 1;
mesh.vtx          = mesh.vtx(logical(Ivtx),:);

% Reorder elements
Ivtx(Ivtx==1) = 1:sum(Ivtx,1);
if size(mesh.elt,1) == 1
    mesh.elt = Ivtx(mesh.elt)';
else
    mesh.elt = Ivtx(mesh.elt);
end
end
