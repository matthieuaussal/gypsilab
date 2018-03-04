%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
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
%|    #    |   FILE       : nrtMshClean.m                                 |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Clean tetrahedral mesh                        |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openMsh')

% Load mesh
load tetmesh

% Build mesh
mesh1 = msh(X,tet);

% Graphical representation
figure
plot(mesh1)
axis equal
view(45,45)

% Add vertices to tetrahedral mesh
tmp = [zeros(10,3) ; 10*rand(10,3)];
tet = tet + size(tmp,1);
X   = [tmp ; X ; 10*rand(20,3)];

% Build and clean mesh
mesh2 = msh(X,tet);

% Test mesh egality
~isequal(mesh1,mesh2)

% Graphical representation
figure
plot(mesh2)
axis equal
view(45,45)



disp('~~> Michto gypsilab !')

