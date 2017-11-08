%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2015-2017.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtFemContinuity.m                            |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2017                                    |
%| ( === ) |   SYNOPSIS   : Continuity beetween domain for a disk         |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openDom')
addpath('../../openFem')
addpath('../../openMsh')

% Parameters
Nvtx = 1e3;
Neig = 10;

% Meshes
mesh = mshDisk(Nvtx,1);

% Colors
ctr = mesh.ctr;
mesh.col(:) = 1;
mesh.col(ctr(:,1)>0) = 2;
mesh.col(ctr(:,1)>0 & ctr(:,2)>0) = 3;

% Submeshing
mesh1 = mesh.sub(mesh.col==1);
mesh2 = mesh.sub(mesh.col==2);
mesh3 = mesh.sub(mesh.col==3);

% Boundary
bound = mesh.bnd;

% Interfaces
inter1 = setdiff(mesh1.bnd,bound);
inter2 = setdiff(mesh2.bnd,bound);

% Domain
omega = dom(mesh,3);
sigma = dom(bound,2);

% Finites elements space
u = fem(mesh,'P1');
v = fem(mesh,'P1');

% Dirichlet
u = dirichlet(u,bound);
v = dirichlet(v,bound);

% Continuity
c = [1 2 ; 1 3 ; 2 3 ];
u = continuity(u,c);
v = continuity(v,c);

% Graphical representation
plot(mesh); 
hold on
plot(bound,'r')
plot(inter1,'y')
plot(inter2,'y')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)

% Mass matrix
tic
M = integral(omega,u,v);
toc

% Rigidity matrix
tic
K = integral(omega,grad(u),grad(v));
toc

% Find eigen values
tic
[V,EV] = eigs(K,M,2*Neig,'SM');
toc

% Evaluate eigen values at vertices
V = feval(u,V,mesh);

% Normalization
V = V./(max(max(abs(V))));

% Sort by ascending order
[EV,ind] = sort(sqrt(real(diag(EV))));
[~,uni]  = unique(floor(1e2*EV));
V        = V(:,ind);

% Graphical representation
figure
for n = 1:9
    subplot(3,3,n)
    plot(mesh,V(:,n))
        title(['k = ',num2str(EV(n))])
    axis equal off
end

% Analytical solutions of eigenvalues for an unitary disk
ref = zeros(Neig^2,1);
k= 1;
for i = 0:Neig
    for j = 0:Neig
        ref(k) = fzero(@(X) besselj(i,X),j);
        k = k+1;
    end
end
ref = ref(ref>1e-6);
ref = 1e-6*unique(ceil(ref*1e6));
ref = ref(1:Neig);

% Error
sol = EV(uni(1:length(ref)));
[ref sol abs(sol(1:length(ref))-ref)./ref]



disp('~~> Michto gypsilab !')


