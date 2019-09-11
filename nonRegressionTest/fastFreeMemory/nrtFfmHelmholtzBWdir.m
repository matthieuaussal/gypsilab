%+========================================================================+
%|                                                                        |
%|         OPENFFM - LIBRARY FOR FAST AND FREE MEMORY CONVOLUTION         |
%|           openFfm is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2019.                             |
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
%|    #    |   FILE       : nrtFfmHelmholtzBWdir.m                        |
%|    #    |   VERSION    : 0.6                                           |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2019                                    |
%| ( === ) |   SYNOPSIS   : Non regression test for BEM resolution of     |
%|  `---'  |                a spherical dirichlet scattering problem      |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parameters
N    = 1e3
tol  = 1e-3
Ninf = 1e2
X0   = [0 0 -1]

% Spherical mesh
mesh = mshSphere(N,1);

% Graphical representation
figure
plot(mesh)
axis equal

% Frequency adjusted to maximum edge size
stp = mesh.stp;
k   = 1/stp(2);
c   = 340;
f   = (k*c)/(2*pi);
disp(['Frequency : ',num2str(f),' Hz']);

% Domain
sigma = dom(mesh,3);    

% Finite elements
u = fem(mesh,'P1');

% Incident wave
PW = @(X) exp(1i*k*X*X0');

% Incident wave representation
hold on
plot(mesh,real(PW(mesh.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
hold off
% view(10,80)
view(0,0)
% camlight
% material dull
% lighting phong


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATORS ~~~~~~~~~~~~~')

% Coupling coeff
beta = 1i*k*0.5;

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(sigma,u,u);

% Finite element boundary regularzation --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
Sr = 1/(4*pi) .* regularize(sigma,sigma,u,'[1/r]',u);
toc

% Finite element boundary regularization --> \int_Sx \int_Sy psi(x)' dny G(x,y) psi(y) dx dy 
tic
Dr = 1/(4*pi) .* regularize(sigma,sigma,u,'grady[1/r]',ntimes(u));
toc

% Close operator : [1i*k*beta*S - (Id/2 + D)]
LHSc = beta.*Sr - (0.5*Id + Dr);

% Left hand side
LHS = @(V) MVdir(sigma,u,k,beta,tol,V) + LHSc*V;

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHS = - integral(sigma,u,PW);


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% LU factorization for preconditionning
tic
[L,U] = ilu(LHSc);
toc

% Solve using iteratve solver
tic
mu     = mgcr(LHS,RHS,[],tol,100,L,U);
lambda = beta * mu;
toc


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy      = '[exp(ikr)/r]';
gradyGxy = {'grady[exp(ikr)/r]1','grady[exp(ikr)/r]2','grady[exp(ikr)/r]3'};

% Mass matrix
Id = integral(sigma,u,u);

% Single layer
tic
Slambda = integral(sigma,sigma,u,Gxy,k,u,tol,lambda);
Slambda = Slambda + regularize(sigma,sigma,u,'[1/r]',u) * lambda;
Slambda = 1/(4*pi) .* Slambda;
toc

% Double layer
tic
Dmu = integral(sigma,sigma,u,gradyGxy,k,ntimes(u),tol,mu);
Dmu = Dmu + regularize(sigma,sigma,u,'grady[1/r]',ntimes(u)) * mu;
Dmu = 1/(4*pi) .* Dmu + 0.5*(Id*mu);
toc

% Boundary solution
Psca = Id \ (Slambda - Dmu);
Pinc = PW(u.dof);
Ptot = Pinc + Psca;

% Graphical representation
figure
plot(mesh,abs(Ptot))
axis equal
title('Total field solution')
colorbar
view(0,0)


%%% ANAYTICAL SOLUTIONS FOR COMPARISONS
% Analytical solution
Ptot = sphereHelmholtz('dom','dir',1,k,1.001*mesh.vtx) + PW(mesh.vtx);

% Solution representation
figure
plot(mesh,abs(Ptot))
axis equal
title('Analytical solution')
colorbar
% view(0,10)


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
theta = 2*pi/Ninf .* (1:Ninf)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty        = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf         = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{1} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,1)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{2} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,2)) .* exp(-1i*k*xdoty(X,Y));
gradxGinf{3} = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,3)) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operators
Sinf = integral(nu,sigma,Ginf,u);
Dinf = integral(nu,sigma,gradxGinf,ntimes(u));

% Finite element radiation  
sol = Sinf*lambda - Dinf*mu;

% Analytical solution
ref = sphereHelmholtz('inf','dir',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')



disp('~~> Michto gypsilab !')


