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
%| Synopsis   : Solve boundary element dirichlet problem, Single Layer    |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openMsh')
addpath('../openHmx')

% Mise en route du calcul paralelle 
% matlabpool; 
% parpool(8)

% Parameters
N   = 1e3
tol = 1e-3
typ = 'P1'
gss = 3
X0  = [0 0 -1]

% Spherical mesh
bnd     = mshSphere(N,1);
bnd.gss = gss;

% Radiative mesh
[x,y] = meshgrid(-3:0.1:3,-3:0.1:3);
DT    = delaunayTriangulation(x(:),y(:));
dom   = msh([zeros(numel(x),1),DT.Points],DT.ConnectivityList);

% Frequency adjusted to maximum esge size
lgt = bnd.edgLgt;
k   = 1/lgt(2)
f   = (k*340)/(2*pi)

% Incident wave
PW = @(X) exp(1i*k*X*X0');

% Incident wave representation
figure
plot(bnd,real(PW(bnd.vtx)))
hold on
plot(dom,real(PW(dom.vtx)))
pts = bnd.ctr;
vct = bnd.nrm;
quiver3(pts(:,1),pts(:,2),pts(:,3),vct(:,1),vct(:,2),vct(:,3),'r');
axis equal;
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
shading interp
hold off


%%% SOLVE LINEAR PROBLEM
disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Finite elements
u = fem(typ);
v = fem(typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
tic
LHS = 1/(4*pi) .* single(integral(bnd,bnd,u,Gxy,v,tol));
toc

% Regularization
tic
LHS = LHS + 1/(4*pi) .* femRegularize(bnd,bnd,u,'[1/r]',v);
toc

% Structure
figure
spy(LHS)

% LU factorization
tic
[Lh,Uh] = lu(LHS);
toc
figure
spy(Lh)
figure
spy(Uh)

% Finite element incident wave trace --> \int_Sx psi(x)' pw(x) dx
RHS = integral(bnd,u,PW);

% Solve linear system [-S] * lambda = - P0
tic
lambda  = Uh \ (Lh \ RHS); % LHS \ RHS;
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
Ninf  = 1e3;
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf  = @(X,Y) 1/(4*pi) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy exp(ik*nu.y) * psi(y) dx
Sinf = integral(nu,bnd,Ginf,v,1e-6);
figure
spy(Sinf)

% Finite element radiation  
sol = - Sinf * lambda;

% Analytical solution
ref = sphereScattering('inf','dir',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
tic
Sdom = 1/(4*pi) .* integral(dom.vtx,bnd,Gxy,v,tol);
toc

% Regularization
tic
Sreg = 1/(4*pi) .* femRegularize(dom.vtx,bnd,'[1/r]',v);
Sdom = Sdom + Sreg;
toc

% Structure
figure
spy(Sdom)

% Boundary solution
Ibnd = integral(bnd,u,v);
Psca = - Ibnd \ double(LHS * lambda) ;
Pinc = PW(u.dof(bnd));
Pbnd = Pinc + Psca;

% Domain solution
Psca = - Sdom * lambda;
Pinc = PW(dom.vtx);
Pdom = Pinc + Psca;

% Annulation sphere interieure
r             = sqrt(sum(dom.vtx.^2,2));
Pdom(r<=1.01) = Pinc(r<=1.01);

% Graphical representation
figure
plot(bnd,abs(Pbnd))
hold on
plot(dom,abs(Pdom))
axis equal;
title('Total field solution')
shading interp
colorbar
hold off


%%% ANAYTICAL SOLUTIONS FOR COMPARISONS
% Analytical solution
Pbnd = sphereScattering('dom','dir',1,k,1.001*bnd.vtx) + PW(bnd.vtx);
Pdom = sphereScattering('dom','dir',1,k,dom.vtx) + PW(dom.vtx);

% Solution representation
figure
plot(bnd,abs(Pbnd))
hold on
plot(dom,abs(Pdom))
axis equal;
title('Analytical solution')
shading interp
colorbar
hold off

disp('Done.')


