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
%| Synopsis   : Eigen values and vectors, disk with dirichlet condition   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openMsh')

% Parameters
Nvtx = 1e3;
Neig = 10;

% Disk mesh
mesh     = mshDisk(Nvtx,1);
mesh.gss = 3;

% Boundary mesh
bnd2vtx = mesh.bnd;
meshb   = msh(mesh.vtx,bnd2vtx);

% Graphical representation
plot(mesh); 
hold on
plot(meshb,'-ob')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)
hold on

% Finites elements
u = fem('P1');
v = fem('P1');

% Mass matrix
tic
M = integral(mesh,u,v);
toc
abs(sum(sum(M)) - pi)/pi

% Rigidity matrix
tic
K = integral(mesh,grad(u),grad(v));
toc
sum(sum(K))

% Dirichlet elimination matrix
P = u.dir(mesh,meshb);

% Eigen value problem matrices
A = P' * K * P;
B = P' * M * P;

% Find eigen values
tic
[V,EV] = eigs(A,B,2*Neig,'SM');
toc

% Sort by ascending order
[EV,ind] = sort(sqrt(real(diag(EV))));
[~,uni]  = unique(floor(1e2*EV));
V        = V(:,ind);

% Eigen vector for all surface and normalisation
V = P*V;
V = V./(max(max(abs(V))));

% Graphical representation
figure
for n = 1:9
    subplot(3,3,n)
    tmp = msh([mesh.vtx(:,1:2) V(:,n)],mesh.elt);
    plot(tmp,V(:,n))
    title(['k = ',num2str(EV(n))])
    shading interp
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


disp('Done, thanks for use.')

