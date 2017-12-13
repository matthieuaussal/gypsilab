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
%|    #    |   FILE       : nrtHmxAlgebra.m                               |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.12.2017                                    |
%| ( === ) |   SYNOPSIS   : Evaluate each hmx class function              |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openHmx')

% Data type
type = 'double';

% Accuracy
tol = 1e-3

% Wave number or frequency (Hz)
k = 1
f = (k*340)/(2*pi);

% Particles receptors X (sphere)
Nx      = 2e3;
[x,y,z] = sphere(ceil(sqrt(Nx)));
X       = unique([x(:),y(:),z(:)],'rows');
Nx      = size(X,1)
if strcmp(type,'single')
    X = single(X);
end

% Particles transmitters Y (=X or not)
Ny = Nx;
Y  = X;
% Ny = 2e3
% Y  = -1+2*rand(Ny,3);
if strcmp(type,'single')
    X = single(X);
end

% Green kernel -> exp(1i*k*r)/r
rxy   = @(X,Y) sqrt( (X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2 + (X(:,3)-Y(:,3)).^2 );
green = @(X,Y) exp(1i*k*rxy(X,Y))./(rxy(X,Y) + 1/(1i*k));

% Particles charges (multiples)
V = (-1+2*rand(Ny,2,type)) + (-1+2i*rand(Ny,2,type));

% Spatial representation of particles
figure
plot3(X(:,1),X(:,2),X(:,3),'*b',Y(:,1),Y(:,2),Y(:,3),'*r')
axis equal 


%%% Full matrix computation
disp('~~~~~~~~~~~~~ EXACT MATRIX AS REFERENCE~~~~~~~~~~~~~')
tic
M = zeros(Nx,Ny,type);
for i = 1:Nx
    M(i,:) = green(X(i,:),Y).';
end
MV = M * V;
toc

tic
un = ones(Nx,1);
I  = spdiags([-un 2*un -un], [-Nx/2,0,Nx/2] , Nx, Ny);
IV = I*double(V);
toc
spy(I)
disp(' ')


%%% H-Matrix computation
disp('~~~~~~~~~~~~~ H-MATRIX ~~~~~~~~~~~~~')
tic
Mh = hmx(X,Y,green,tol);
toc

tic
Mh2 = hmx(X,Y,M,tol);
toc

tic
Ih = hmx(double(X),double(Y),I,tol);
toc

disp(' ')


%%% H-Matrix structure
disp('~~~~~~~~~~~~~ SPY STRUCTURE ~~~~~~~~~~~~~')
tic
figure
spy(Mh);
toc

tic
figure
spy(Mh2);
toc

tic
figure
spy(Ih);
toc

disp(' ')


%%% Full conversion
disp('~~~~~~~~~~~~~ FULL CONVERSION ~~~~~~~~~~~~~')
tic
sol = full(Mh);
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = full(Mh2);
toc
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = full(Ih);
toc
ref = full(I);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Recompression
disp('~~~~~~~~~~~~~ RECOMPRESSION ~~~~~~~~~~~~~')
tic
tmp = hmxRecompress(Mh,10*tol);
toc
sol = full(tmp);
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Sparse conversion
disp('~~~~~~~~~~~~~ SPARSE CONVERSION ~~~~~~~~~~~~~')
tic
sol = sparse(double(Mh));
toc
ref = double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = sparse(Ih);
toc
ref = I;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Low-rank conversion
disp('~~~~~~~~~~~~~ LOW-RANK CONVERSION ~~~~~~~~~~~~~')
tic
[A,B] = hmxLowrank(double(Mh));
sol   = A*B;
toc
ref = double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Concatenation
disp('~~~~~~~~~~~~~ CONCATENATION ~~~~~~~~~~~~~')
AB = hmx(X,Y,A,B,tol);
tic
tmp = [AB,Mh,Mh,AB];
toc
sol = full(tmp);
ref = [M,M,M,M];
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = [AB;Mh;Mh;AB];
toc
sol = full(tmp);
ref = [M;M;M;M];
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = [Ih,Ih];
toc
sol = full(tmp);
ref = [I,I];
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Transposition
disp('~~~~~~~~~~~~~ TRANSPOSITION ~~~~~~~~~~~~~')
tic
tmp = Mh.';
toc
sol = full(tmp);
ref = M.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih.';
toc
sol = full(tmp);
ref = I.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Conjugate transposition
disp('~~~~~~~~~~~~~ CTRANSPOSITION ~~~~~~~~~~~~~')
tic
tmp = Mh';
toc
sol = full(tmp);
ref = M';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih';
toc
sol = full(tmp);
ref = I';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix right product
disp('~~~~~~~~~~~~~ H-MATRIX * MATRIX ~~~~~~~~~~~~~')
tic
sol = Mh * V;
toc
ref = MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = Ih * double(V);
toc
ref = IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = double(Mh) * I.';
toc
ref = double(M) * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol = Ih * I.';
toc
ref = I * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Left product H-Matrix
disp('~~~~~~~~~~~~~ FULL MATRIX * H-MATRIX ~~~~~~~~~~~~~')
tmp = Mh.';
tic
sol = V.' * tmp;
toc
ref = MV.';
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Scalar product
disp('~~~~~~~~~~~~~ SCALAR PRODUCT ~~~~~~~~~~~~~')
tic
tmp = (sqrt(2) .* Mh .* pi);
toc
sol = tmp * V;
ref = sqrt(2) * pi * MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = (sqrt(2) .* Ih .* pi);
toc
sol = tmp * double(V);
ref = sqrt(2) * pi * IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Uminus
disp('~~~~~~~~~~~~~ UMINUS ~~~~~~~~~~~~~')
tic
tmp = - Mh;
toc
sol = tmp * V;
ref = - MV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Full Addition
disp('~~~~~~~~~~~~~ FULL ADDITION/SUBSTRACTION ~~~~~~~~~~~~~')
tic
tmp = 2*M + Mh - M;
toc
sol = tmp * V;
ref = 2*MV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2*M + Ih - M;
toc
sol = tmp * V;
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Full Addition
disp('~~~~~~~~~~~~~ SPARSE ADDITION/SUBSTRACTION ~~~~~~~~~~~~~')
tic
tmp = 2*I + Mh - I;
toc
sol = tmp * V;
ref = MV+IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2*I + Ih - I;
toc
sol = tmp * double(V);
ref = 2*IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix addition
disp('~~~~~~~~~~~~~ H-MATRIX ADDITION/SUBSTRACTION ~~~~~~~~~~~~~')
tic
tmp = 2.*Mh + Mh - Mh;
toc
sol = tmp * V;
ref = 2*MV;
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

tic
tmp = 2.*Ih + Mh - Ih;
toc
sol = tmp * V;
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*Mh + Ih - Mh;
toc
sol = tmp * V;
ref = MV + IV;
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = 2.*Ih + Ih - Ih;
toc
sol = tmp * double(V);
ref = 2*IV;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix product
disp('~~~~~~~~~~~~~ H-MATRIX PRODUCT ~~~~~~~~~~~~~')
tic
tmp = Mh * Mh.';
toc
sol = full(tmp);
ref = M * M.';
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

tic
tmp = Ih.' * double(Mh);
toc
sol = full(tmp);
ref = I.' * double(M);
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = double(Mh) * Ih.';
toc
sol = full(tmp);
ref = double(M) * I.';
norm(ref-sol,'inf')/norm(ref,'inf')

tic
tmp = Ih * Ih.';
toc
sol = sparse(tmp);
ref = I * I.';
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

disp(' ')


%%% H-Matrix inversion
disp('~~~~~~~~~~~~~ INVERSION ~~~~~~~~~~~~~')
Mh = Mh + sqrt(Nx).*speye(Nx);
M  = M  + sqrt(Nx).*eye(Nx,type); 

tic
tmp = inv(Mh);
toc
tic
ref = inv(M);
toc
sol = full(tmp);
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

tic
tmp = inv(Ih);
toc
tic
ref = inv(I);
toc
sol = full(tmp);
norm(ref-sol,'inf')/norm(ref,'inf')
figure
spy(tmp)

disp(' ')


%%% H-Matrix cholesky
disp('~~~~~~~~~~~~~ CHOLESKY FACTORISATION ~~~~~~~~~~~~~')
tic
Uh = chol(Ih);
toc
sol = sparse(Uh'*Uh);
tic
U = chol(I);
toc
sol2 = U'*U;
ref  = I;
norm(ref-sol,'inf')/norm(ref,'inf')
norm(ref-sol2,'inf')/norm(ref,'inf')
figure
spy(Uh)

disp(' ')


%%% H-Matrix LDLt
disp('~~~~~~~~~~~~~ LDLt FACTORISATION ~~~~~~~~~~~~~')
tic
[Lh,Dh] = ldl(Ih);
toc
sol = sparse(Lh*Dh*Lh');
tic
[L,D,P] = ldl(I);
L       = P*L;
toc
sol2 = L*D*L';
ref  = I;
norm(ref-sol,'inf')/norm(ref,'inf')
norm(ref-sol2,'inf')/norm(ref,'inf')
figure
spy(Lh)

disp(' ')


%%% H-Matrix LU
disp('~~~~~~~~~~~~~ LU FACTORISATION ~~~~~~~~~~~~~')
tic
[Lh,Uh] = lu(Mh);
toc
tmp = full(Mh);
tic
[L,U] = lu(tmp);
toc
sol  = full(Lh*Uh);
sol2 = L*U;
ref  = M;
norm(ref-sol,'inf')/norm(ref,'inf')
norm(ref-sol2,'inf')/norm(ref,'inf')
figure
spy(Lh)
figure
spy(Uh)

tic
[Lh,Uh] = lu(Ih);
toc
sol = full(Lh*Uh);
tic
[L,U] = lu(I);
toc
sol2 = L*U;
ref  = I;
norm(ref-sol,'inf')/norm(ref,'inf')
norm(ref-sol2,'inf')/norm(ref,'inf')
figure
spy(Lh)
figure
spy(Uh)

disp(' ')


%%% Shermann-morrisonn conversion
disp('~~~~~~~~~~~~~ SHERMANN MORRISON CONVERSION ~~~~~~~~~~~~~')
tic
[Sh,A,B] = hmxSherMorr(Mh);
toc
sol = full(Sh) + A*B;
ref = M;
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% H-Matrix \
disp('~~~~~~~~~~~~~ LU SOLVER ~~~~~~~~~~~~~')
tic
sol =  Mh \ V;
toc
tic
ref = M \ V;
toc
norm(ref-sol,'inf')/norm(ref,'inf')

tic
sol =  Ih \ double(V);
toc
tic
ref = I \ double(V);
toc
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% Shermann-morrisonn solver
disp('~~~~~~~~~~~~~ SHERMANN MORRISON SOLVER ~~~~~~~~~~~~~')

tic
[Lh,Uh] = lu(Sh);
toc

tic
rk   = size(A,2);
Sm1  = Uh \ (Lh \ [A,V]);
Sm1A = Sm1(:,1:rk);
Sm1V = Sm1(:,rk+1:end);
toc

tic
Mk    = eye(rk) + B*(Sm1A);
Mkm1V = Mk \ (B * Sm1V);
toc

tic
sol = Sm1V - Uh \ (Lh \ ( A * Mkm1V));
toc
ref = M \ V;

norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% ITERATIVE SOLVER
disp('~~~~~~~~~~~~~ ITERATIVE SOLVER ~~~~~~~~~~~~~')
tic
sol = zeros(Nx,size(V,2),type);
for i = 1:size(V,2)
    sol(:,i) = gmres(@(V) Mh*V,V(:,i),[],tol,100);
end
toc
ref = zeros(Nx,size(V,2),type);
for i = 1:size(V,2)
    ref(:,i) = gmres(M,V(:,i),[],tol,100);
end
toc
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')


%%% ITERATIVE SOLVER WITH PRECONDITIONNER
disp('~~~~~~~~~~~~~ ITERATIVE SOLVER WITH PRECONDITIONNER ~~~~~~~~~~~~~')
tic
[Lh,Uh] = lu(hmxRecompress(Mh,0.1));
Mhm1V   = @(V) Uh\(Lh\V);
toc
tic
sol = zeros(Nx,size(V,2),type);
for i = 1:size(V,2)
    sol(:,i) = gmres(@(V) Mh*V,V(:,i),[],tol,100,Mhm1V);
end
toc
norm(ref-sol,'inf')/norm(ref,'inf')

disp(' ')




disp('~~> Michto gypsilab !')







% % Preconditionner
% S  = sparse(Sh);
% Ap = A(:,1:1);
% Bp = B(1:1,:);
% D  = eye(1) + Bp*(S\Ap);
% 
% P = @(X) S\X - S\(Ap*(D\(Bp*(S\X))));
% 
% norm(S*V + Ap*(Bp*V) - MV,'inf')/norm(MV,'inf')
% norm(P(V) - ref,'inf')/norm(ref,'inf')
% 
% tic
% sol = zeros(Nx,size(V,2),type);
% for i = 1:size(V,2)
%     sol(:,i) = gmres(@(V) Mh*V,V(:,i),[],1e-3,100,@(X));
% end
% toc
% ref = zeros(Nx,size(V,2),type);
% for i = 1:size(V,2)
%     ref(:,i) = gmres(M,V(:,i),[],tol,100);
% end
% toc
% norm(ref-sol,'inf')/norm(ref,'inf')

