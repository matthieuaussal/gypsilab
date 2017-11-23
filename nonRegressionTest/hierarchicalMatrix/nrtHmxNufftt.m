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
%|    #    |   FILE       : nrtHmxBuilder.m                               |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Non uniform FFT with ACA and H-Matrix         |
%|  `---'  |                                                              |
%+========================================================================+

% Nettoyage
clear all
close all
clc

% Library path
addpath('../../openFem')
addpath('../../openHmx')

% Mise en route du calcul paralelle 
% matlabpool; 
% parpool

% Spheres unitaires
unit = @(u) u ./ sqrt(sum(u.^2,2)*ones(1,size(u,2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Type de donnees
type = 'double'

% Precision
tol = 1e-3

% Nombre d'onde et frequence (Hz)
k = 10
f = (k*340)/(2*pi);

% Nuage de point recepteurs X
Nx = 1e4
X  = -1 + 2*rand(Nx,3,type);
X  = unit(X);

% Nuage de point emmeteurs Y
Ny = 1e4
Y  = -1 + 2*rand(Nx,3,type);
% Y  = unit(Y);

% Potentiel en Y
V = -(1+1i) + 2*(rand(Ny,1,type) + 1i*rand(Ny,1,type));

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
xdoty = @(X,Y) X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3);
Gxy   = @(X,Y) exp(1i*k*xdoty(X,Y));
% rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
% Gxy   = @(X,Y) exp(1i*k*rxy(X,Y))./(rxy(X,Y) + 1/(1i*k));

% Representation graphique
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCUL DIRECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ FULL PRODUCT ~~~~~~~~~~~~~')
tic

% Recepteurs au hasard
Nt   = 100;
ind  = 1:ceil(Nx/Nt):Nx;
Xloc = X(ind,:);

% Produit Matrice-Vecteur sur tous les emeteurs
ref = zeros(length(ind),1,type);
for i = 1:length(ind)
    ref(i) = Gxy(Xloc(i,:),Y).' * V;
end
toc
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPRESSION ACA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ ACA COMPRESSION ~~~~~~~~~~~~~')
tic

% Compression ACA
tic
[A,B] = hmxACA(X,Y,Gxy,tol);
toc
size(A)

% Erreur
tic
sol = A * (B * V);
toc
norm(ref-sol(ind))/norm(ref)
disp(' ')


%%%%%%%%%%%%%%%%%%%%%% COMPRESSION ACA HIERARCHIQUE %%%%%%%%%%%%%%%%%%%%%%%
disp('~~~~~~~~~~~~~ H-MATRIX COMPRESSION ~~~~~~~~~~~~~')

% Hierarchical Matrix
tic
Mh = hmx(X,Y,Gxy,tol);
toc

% Size
% tmp = whos('Mh');
% tmp.bytes/1e6

% Graphical representation
figure
spy(Mh)

% Matrix vector product
tic
sol = Mh * V;
toc

% Erreur
norm(ref-sol(ind))/norm(ref)
disp(' ')


disp('~~> Michto gypsilab !')

