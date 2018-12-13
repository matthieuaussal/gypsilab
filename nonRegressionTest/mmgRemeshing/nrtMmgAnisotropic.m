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
%|    #    |   FILE       : nrtMmgAnisotropic.m                           |
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

% Homothetie
mesh.vtx = 10+100*mesh.vtx;

% Define anisotropic refinment
X   = mesh.vtx;
stp = mesh.stp
in  = stp(3) * ones(size(X,1),1);
in(abs(X(:,2))<=1) = 1;

% Mesh refinment
[meshMmg,out] = mmg(mesh,in);
meshMmg.stp
    
% Graphical representation
figure

subplot(1,2,1)

plot(mesh,'w')
hold on
plot(mesh,in)
alpha(0.5)
axis equal
view(0,90)    
grid on
title('Original')
colorbar

subplot(1,2,2)
plot(meshMmg,'w')
hold on
plot(meshMmg,out)
alpha(0.5)
axis equal
view(0,90)    
grid on
title('Mmg')
colorbar


%%% 3D surface
N = 1000;

% Mesh unit sphere
mesh = mshSphere(N,1);

% Homothetie
mesh.vtx = 10+100*mesh.vtx;

% Define anisotropic refinment
X   = mesh.vtx;
stp = mesh.stp
in  = stp(3) * ones(size(X,1),1);
in(abs(X(:,2))<=10) = 4;

% Mesh refinment
[meshMmg,out] = mmg(mesh,in);
meshMmg.stp
    
% Graphical representation
figure

subplot(1,2,1)
plot(mesh)
hold on
plot(mesh,in,'w')
axis equal
view(0,90)    
grid on
title('Original')
colorbar

subplot(1,2,2)
plot(meshMmg,'w')
hold on
plot(meshMmg,out)
axis equal
view(0,90)    
grid on
title('Mmg')
colorbar


%%% 3D volume
N = 100;

% Mesh unit cube
mesh = mshCube(N,[1 1 1]);

% Homothetie
mesh.vtx = 10+100*mesh.vtx;

% Define anisotropic refinment
X   = mesh.vtx;
stp = mesh.stp
in  = stp(3) * ones(size(X,1),1);
in(abs(X(:,3))<=5) = 5;

% Mesh refinment
[meshMmg,out] = mmg(mesh,in);
meshMmg.stp

% Graphical representation
figure

subplot(1,2,1)
plot(mesh.edg,'k')
hold on
plot(mesh,in)
alpha(0.5)
axis equal
view(10,10)    
grid on
title('Original')
colorbar

subplot(1,2,2)
plot(meshMmg.edg,'k')
hold on
plot(meshMmg,out)
alpha(0.5)  
axis equal
view(10,10)    
grid on
title('Mmg')
colorbar

    

disp('~~> Michto gypsilab !')



