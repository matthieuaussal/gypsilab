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
%| Synopsis   : Solve boundary element dirichlet problem, Double Layer    |
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
k   = 1/lgt(2);
f   = (k*340)/(2*pi);

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

% Finite elements
u = fem(typ);
v = fem(typ);

% Finite element mass matrix --> \int_Sx psi(x)' psi(x) dx
Id = integral(bnd,u,v);

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy1 = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k);
Gxy2 = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k);
Gxy3 = @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' grady(G(x,y)) ny.psi(y) dx dy 
tic
Dbnd = 1/(4*pi) .* (integral(bnd,bnd,u,Gxy1,ndot(v,1),tol) + ...
    integral(bnd,bnd,u,Gxy2,ndot(v,2),tol) + ...
    integral(bnd,bnd,u,Gxy3,ndot(v,3),tol) ) ;
toc

% Regularization
tic
Dbnd = Dbnd + 1/(4*pi) .* (femRegularize(bnd,bnd,u,'grady[1/r]1',ndot(v,1)) + ...
    femRegularize(bnd,bnd,u,'grady[1/r]2',ndot(v,2)) + ...
    femRegularize(bnd,bnd,u,'grady[1/r]3',ndot(v,3)) ) ;
toc

% Operator [Id/2 + D]
LHS = 0.5*Id + Dbnd;

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

% Finite element incident wave trace --> \int_Sx psi(x) pw(x) dx
RHS = - integral(bnd,u,PW);

% Solve linear system [Id/2 + D] * mu = - P0
tic
mu  = Uh \ (Lh \ RHS); % LHS \ RHS;
toc


%%% INFINITE SOLUTION
disp('~~~~~~~~~~~~~ INFINITE RADIATION ~~~~~~~~~~~~~')

% Plane waves direction
Ninf  = 1e3;
theta = 2*pi/1e3 .* (1:1e3)';
nu    = [sin(theta),zeros(size(theta)),cos(theta)];

% Green kernel function
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3); 
Ginf1 = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,1)) .* exp(-1i*k*xdoty(X,Y));
Ginf2 = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,2)) .* exp(-1i*k*xdoty(X,Y));
Ginf3 = @(X,Y) 1/(4*pi) .* (-1i*k*X(:,3)) .* exp(-1i*k*xdoty(X,Y));

% Finite element infinite operator --> \int_Sy dny(exp(ik*nu.y)) * psi(y) dx
Dinf = integral(nu,bnd,Ginf1,ndot(v,1),1e-6) + ...
    integral(nu,bnd,Ginf2,ndot(v,2),1e-6) + ...
    integral(nu,bnd,Ginf3,ndot(v,3),1e-6) ;
figure
spy(Dinf)

% Finite element radiation  
sol = Dinf * mu;

% Analytical solution
ref = sphereScattering('inf','dir',1,k,nu); 
norm(ref-sol,2)/norm(ref,2)
norm(ref-sol,'inf')/norm(ref,'inf')

% Graphical representation
figure
plot(theta,log(abs(sol)),'b',theta,log(abs(ref)),'--r')


%%% DOMAIN SOLUTION
disp('~~~~~~~~~~~~~ RADIATION ~~~~~~~~~~~~~')

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' grady(G(x,y)) ny.psi(y) dx dy 
tic
Ddom = 1/(4*pi) .* (integral(dom.vtx,bnd,Gxy1,ndot(v,1),tol) + ...
    integral(dom.vtx,bnd,Gxy2,ndot(v,2),tol) + ...
    integral(dom.vtx,bnd,Gxy3,ndot(v,3),tol) ) ;
toc

% Regularization
tic
Ddom = Ddom + 1/(4*pi) .* (femRegularize(dom.vtx,bnd,'grady[1/r]1',ndot(v,1)) + ...
    femRegularize(dom.vtx,bnd,'grady[1/r]2',ndot(v,2)) + ...
    femRegularize(dom.vtx,bnd,'grady[1/r]3',ndot(v,3)) ) ;
toc

% Structure
figure
spy(Ddom)

% Boundary solution
Psca = 0.5*mu + Id \ (Dbnd * mu) ;
Pinc = PW(u.dof(bnd));
Pbnd = Pinc + Psca;

% Domain solution
Psca = Ddom * mu;
Pinc = PW(dom.vtx);
Pdom = Pinc + Psca;

% Annulation sphere interieure
r              = sqrt(sum(dom.vtx.^2,2));
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
