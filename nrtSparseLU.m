%|                                                                        |
%|                      THE GYPSILAB TOOLBOX - v0.20                      |
%|               www.cmap.polytechnique.fr/~aussal/gypsilab               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. You can         \ \| |/ /   |
%| refer to each library part each part for references.        \ | | /    |
%|                                                              \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Comparse sparse LU with H-Matrix sparse LU                |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openHmx')
addpath('../openMsh')

% Parameters
Nvtx = 1e4;
L    = [1 0.1 0.1];

% Cube mesh
mesh     = mshCube(Nvtx,L);
mesh.gss = 4;

% Graphical representation
plot(mesh); 
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.1)
hold on

% Finites elements
u = fem('P1');
v = fem('P1');

% Random vector at dof
V = rand(size(v.dof(mesh),1),1);

% Mass matrix
tic
M = integral(mesh,u,v);
toc
ref = M * V;
figure
spy(M)
drawnow

% LU factorisation
tic
[L,U] = lu(M);
toc
norm(L*(U*V) - ref)/norm(ref)
figure
spy(L)

% H-Matrix
tic
Mh = hmx(u.dof(mesh),v.dof(mesh),M,1e-12);
toc
norm(Mh*V - ref)/norm(ref)
figure
spy(Mh)

% LhUh factorisation
tic
[Lh,Uh] = lu(Mh);
toc
norm(Lh*(Uh*V) - ref)/norm(ref)
figure
spy(Lh)
drawnow

% Solver
tic
Mm1V = U \ (L \ V);
toc
tic
Mhm1V = Uh \ (Lh \ V);
toc
tic
ref = M \ V;
toc
norm(Mm1V - ref)/norm(ref)
norm(Mhm1V - ref)/norm(ref)

whos


