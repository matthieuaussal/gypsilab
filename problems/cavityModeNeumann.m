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
%| Synopsis   : Exterior eigen mode for a cubic cavity, neumann condition |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openHmx')
addpath('../openMsh')

% Mise en route du calcul paralelle 
% matlabpool; 
% parpool(8)

% Parameters
X0  = [0 0 -1];
f   = 300;
typ = 'P1';
tol = 1e-3;

% Spherical mesh
bnd     = msh('meshes/unitCubeCavity_1e3.ply');
bnd.gss = 3;

% Radiating mesh
dom = msh('meshes/Oxz_1e4.ply');

% Wave number with high frequency security
k = f*2*pi/340;

% Incident wave
PW = @(X) exp(1i*k*X*X0');

% Incident wave representation
figure
plot(bnd,real(PW(bnd.vtx)))
alpha(0.7)
hold on
plot(dom,real(PW(dom.vtx)))
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

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Finite element boundary operator --> 
% k^2 * \int_Sx \int_Sy n.psi(x) G(x,y) n.psi(y) dx dy 
% - \int_Sx \int_Sy nxgrad(psi(x)) G(x,y) nxgrad(psi(y)) dx dy 
tic
LHS = 1/(4*pi) .* (k^2 * single(integral(bnd,bnd,ndot(u),Gxy,ndot(v),tol)) ...
    - single(integral(bnd,bnd,nxgrad(u),Gxy,nxgrad(v),tol)));
toc

% Regularization
tic
LHS = LHS + 1/(4*pi) .* (k^2 * femRegularize(bnd,bnd,ndot(u),'[1/r]',ndot(v)) ...
    - femRegularize(bnd,bnd,nxgrad(u),'[1/r]',nxgrad(v)));
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

% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
gradxPW1 = @(X) 1i*k*X0(1) .* PW(X);
gradxPW2 = @(X) 1i*k*X0(2) .* PW(X);
gradxPW3 = @(X) 1i*k*X0(3) .* PW(X);
RHS = - (integral(bnd,ndot(u,1),gradxPW1) + ...
    integral(bnd,ndot(u,2),gradxPW2) + ...
    integral(bnd,ndot(u,3),gradxPW3) );

% Solve linear system  [H] mu = - dnP0
tic
mu = Uh \ (Lh \ RHS); % LHS \ RHS;
toc


%%% BOUNDARY RADIATION 
disp('~~~~~~~~~~~~~ RADIATION ON THE BOUNDARY ~~~~~~~~~~~~~')

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

% Boundary solution
Psca = 0.5*mu + Id \ double(Dbnd * mu) ;
Pinc = PW(u.dof(bnd));
Pbnd = Pinc + Psca;

% Graphical representation
figure
plot(bnd,abs(Pbnd))
alpha(0.7)
axis equal;
title('Total field solution')
shading interp
colorbar


%%% DOMAIN RADIATION 
disp('~~~~~~~~~~~~~ RADIATION ON THE DOMAIN ~~~~~~~~~~~~~')

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

% Boundary solution
Psca = double(Ddom * mu) ;
Pinc = PW(dom.vtx);
Pdom = Pinc + Psca;

% Graphical representation
figure
plot(dom,abs(Pdom))
axis equal;
title('Total field solution')
shading interp
colorbar


disp('Done.')
