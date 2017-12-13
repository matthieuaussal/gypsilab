%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2015-2017.                             |
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
%|    #    |   FILE       : nrtHmxBuilderCIVA.m                           |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Build H-Matrix for CIVA                       |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Library path
addpath('../../openMsh')
addpath('../../openDom')
addpath('../../openFem')
addpath('../../openHmx')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametres
N   = 1e3;
fe  = 'P1';
gss = 3;
tol = 1e-3;

% Mesh
mesh = mshSphere(N,1);

% Nombre d'onde et frequence (Hz)
stp = mesh.stp;
k   = 0;%1/stp(2);
f   = (k*340)/(2*pi);

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',k);

% Quadrature
omega = dom(mesh,gss);
X     = omega.qud;
Y     = X;

% Elements finis
u    = fem(mesh,fe);
DQx  = u.dqm(omega)';
QDy  = DQx';

% Inconnues
ctr  = mesh.ctr;
half = mesh.sub(ctr(:,1)>=-0.5);
UDx  = restriction(u,half);
DUy  = UDx';

% Interactions boolennes inconnues <-> quadratures 
Bx = (abs(UDx*DQx) > 0);
By = (abs(QDy*DUy) > 0)';

% Interactions proches
Dxy = integral(omega,u,u);

% Potentiel en Uy
V = -(1+1i) + 2*(rand(size(UDx,1),1) + 1i*rand(size(UDx,1),1));

% Graphical representation
figure
plot(mesh)
hold on
plot(half,'r')
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')

% Full matrix at gauss points
tic
Mxy = zeros(size(X,1),size(Y,1));
for j = 1:size(Y,1)
    Mxy(:,j) = Gxy(X,Y(j,:));
end
toc

% Quadrature -> dof
M = DQx * Mxy * QDy;

% Replace close operator
[I,J,S] = find(Dxy);
M(sub2ind(size(M),I,J)) = S;

% Dof -> unknowns
M = UDx * M * DUy;

% Matrix-vector product
ref = M*V;

% Facorisation LU
tic
[L,U] = lu(M);
toc
sol = L*(U*V);

% Erreur
norm(ref-sol)/norm(ref)
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%% ORIGINAL CONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ ORIGINAL CONSTRUCTION ~~~~~~~~~~~~~')

'FULL'

% Hierarchical Matrix
tic
Dh = hmx(u.dof,u.dof,DQx,X,Gxy,Y,QDy,tol); %%%%%% 
toc

% Correct close part
tic
Dh = hmxInsert(Dh,Dxy);
toc

% Graphical representation
figure
spy(Dh)

'SPARSE'

% Hierarchical sparse matrix
tic
Uh = hmxBuilder(UDx*u.dof,u.dof,UDx,tol);
toc

% Graphical representation
figure
spy(Uh)

'PRODUCT'

% Dof -> unknowns
tic
Mh = Uh * Dh;
toc
tic
Mh = Mh * Uh.';
toc

% Graphical representation
figure
spy(Mh)

% Matrix vector product
tic
sol = Mh * V;
toc

% Erreur
norm(ref-sol)/norm(ref)
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LU FACTORIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ LU FACTORIZATION ~~~~~~~~~~~~~')

% Factorization LU
tic
[Lh,Uh] = lu(Mh);
toc

% Graphical representation
figure
subplot(1,2,1)
spy(Lh)
subplot(1,2,2)
spy(Uh)

% Matrix vector product
sol = Lh * (Uh * V);
norm(ref-sol)/norm(ref)
disp(' ')


disp('~~> Michto gypsilab !')

