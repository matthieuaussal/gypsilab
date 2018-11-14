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
%|    #    |   FILE       : nrtStkHmx.m                                   |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2018                                    |
%| ( === ) |   SYNOPSIS   : Stoke Galerkin of a unit sphere using H-matrix|
%|  `---'  |                with stokeslet G and stresslet T              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openDom')
addpath('../../openFem')
addpath('../../openMsh')
addpath('../../openHmx')

% Parameters
N   = 100;
x0  = [0.1 0.2 0.3];
tol = 1e-3;

% Mesh unit sphere
mesh = mshSphere(N,1);

% Quadrature
gamma = dom(mesh,3);

% Finite element
phi = fem(mesh,'P1');
unk = phi.unk;
un  = ones(length(phi),1);

% Mass matrix
M = integral(gamma,phi,phi);

% Graphical rep
figure
plot(mesh)
hold on
plotNrm(mesh)
plot(gamma)
plot(phi,'*r')
axis equal

% Initialization with cells
mu     = cell(3,1);
lambda = cell(3,1);
G      = cell(6,1);
T      = cell(6,1);
C      = cell(6,1);
I      = cell(6,1);
ind    = { [1 1] ; [2 2] ; [3 3] ; [1 2] ; [1 3] ; [2 3] };

% Vectors
tic
for i = 1:3
    % mu = [u] = u_int - u_ext = - Gi1 = - Gi1(x0,y)
    name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(1)];
    green = @(Y) 1/(8*pi) .* femGreenKernel(x0,Y,name,[]);
    mu{i} = - integral(gamma,phi,green);

    % lambda = [sigma] = - Ti1
    lambda{i} = 0;
    for k = 1:3
        name      = ['[rirjrk/r^5]',num2str(i),num2str(1),num2str(k)];
        green     = @(Y) -6/(8*pi) .* femGreenKernel(x0,Y,name,[]);
        lambda{i} = lambda{i} - integral(gamma,ntimes(phi,k),green);
    end
end
toc

% Matrices using symetry
tic
for n = 1:6
    % Indices
    i = ind{n}(1);
    j = ind{n}(2);
    
    % Single layer : G = \int_gamma \int_gamma  1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3)
    name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
    green = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
    G{n}  = integral(gamma,gamma,phi,green,phi,tol);
    
    % Regularization
    G{n} = G{n} + 1/(8*pi) .* regularize(gamma,gamma,phi,name,phi);
    
    % Double layer  : T = \int_gamma \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
    % Double lumped : P = \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
    T{n} = zeros(hmx(unk,unk,tol));
    P    = 0;
    for k = 1:3
        name   = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
        green  = @(X,Y) -6/(8*pi) .* femGreenKernel(X,Y,name,[]);
        T{n}   = T{n} + integral(gamma,gamma,phi,green,ntimes(phi,k),tol);
        P      = P + integral(gamma.qud,gamma,green,ntimes(phi,k),tol) * un;
    end
    
    % Correction
    C{n} = integral(gamma,phi,@(X)P,phi);
    
    % Identity matrix
    if (i == j)
        I{n} = M;
    else
        I{n} = sparse(size(M,1),size(M,2));
    end
end
toc

% Convert cells to matrix
tic
G = [G{1} , G{4} , G{5} ; G{4} , G{2} , G{6} ; G{5} , G{6} , G{3} ];
T = [T{1} , T{4} , T{5} ; T{4} , T{2} , T{6} ; T{5} , T{6} , T{3} ];
C = [C{1} , C{4} , C{5} ; C{4} , C{2} , C{6} ; C{5} , C{6} , C{3} ];
I = [I{1} , I{4} , I{5} ; I{4} , I{2} , I{6} ; I{5} , I{6} , I{3} ];
toc
mu     = cell2mat(mu);
lambda = cell2mat(lambda);

% Graphical representation
figure
subplot(1,2,1)
spy(G)
subplot(1,2,2)
spy(T)

% Validation test (cf. Aline Lefebvre-Lepot)
sum(mu) - (-2/3)
sum(lambda) - (-1)

% Analytic solution : Gi1
ref = - mu;

% Stokes radiation in boundary : ui(x) = - sum_j \int_gamma Gij(x,y) lambda_j dy ...
%              + sum_j \int_gamma Tij(x,y).n(y) mu_j dy 
sol = -G*(I\lambda) + T*(I\mu) - 0.5*mu; 
norm(ref-sol)/norm(ref)
norm(ref-sol,'inf')/norm(ref,'inf')

% Stokes radiation with correction
sol = -G*(I\lambda) + T*(I\mu) - C*(I\mu); 
norm(ref-sol)/norm(ref)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation of solution
figure
for i = 1:3
    subplot(1,3,i)
    ind = (i-1)*N + (1:N);
    plot(mesh,sol(ind))
    axis equal
    colorbar
    title(['Component ',num2str(i)])
end





disp('~~> Michto gypsilab !')




