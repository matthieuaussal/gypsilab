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
%|    #    |   FILE       : nrtMmgOptim.m                                 |
%|    #    |   VERSION    : 0.50                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 25.11.2018                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   : Remeshing mesh using Mmg platform             |
%|  `---'  |                https://www.mmgtools.org/                     |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openMsh')
addpath('../../openMmg')


%%% 2D surface
N = 100;

% Mesh unit square
mesh = mshSquare(N,[1 1]);

% Move points
mesh.vtx(:,1) = 10*mesh.vtx(:,1); 
mesh.stp

% Mesh refinment
meshMmg = mmg(mesh);
meshMmg.stp
    
% Graphical representation
figure

subplot(1,2,1)
plot(mesh)
axis equal
view(0,90)    
grid on
title('Original')

subplot(1,2,2)
plot(meshMmg)
axis equal
view(0,90)    
grid on
title('Mmg')


%%% 3D surface
N = 100;

% Mesh unit sphere
mesh = mshSphere(N,1);

% Move points
mesh.vtx(:,1) = 5*mesh.vtx(:,1); 
mesh.stp

% Mesh refinment
meshMmg = mmg(mesh);
meshMmg.stp
    
% Graphical representation
figure

subplot(1,2,1)
plot(mesh)
axis equal
view(0,90)    
grid on
title('Original')

subplot(1,2,2)
plot(meshMmg)
axis equal
view(0,90)    
grid on
title('Mmg')


%%% 3D volume
N    = 100;

% Mesh unit cube
mesh = mshCube(N,[1 1 1]);

% Move points
mesh.vtx(:,1) = 10*mesh.vtx(:,1); 
mesh.stp

% Mesh refinment
meshMmg = mmg(mesh);
meshMmg.stp

% Graphical representation
figure

subplot(1,2,1)
plot(mesh)
axis equal
view(20,20)    
grid on
title('Original')
alpha(0.3)

subplot(1,2,2)
plot(meshMmg)
axis equal
view(20,20)    
grid on
title('Mmg')
alpha(0.3)
    

    

disp('~~> Michto gypsilab !')


