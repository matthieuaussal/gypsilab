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
%| Synopsis   : Eigen values and vectors, cube with dirichlet condition   |
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
L    = [1 0.5 0.5];

% Cube mesh
mesh = mshCube(Nvtx,L);
% load tetmesh
% mesh = msh(X,tet);
mesh.gss = 4;

% Boundary mesh
bnd2vtx = mesh.bnd;
meshb   = msh(mesh.vtx,bnd2vtx);

% Half boundary
ctr   = mesh.ctr;
Iplus = find(ctr(:,1)>0);
meshh = mesh.sub(Iplus);

% Graphical representation
plot(mesh); 
hold on
plot(meshh,'r')
hold off
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.1)
hold on

% Finites elements
u = fem('P1');
v = fem('P1');

% Mass matrix
tic
M = integral(mesh,u,v);
toc
abs(sum(sum(M)) - prod(L))/prod(L)

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
V        = V(:,ind);

% Eigen vector for all surface and normalisation
V = P*V;
V = V./(max(max(abs(V))));

% Graphical representation
figure
for n = 1:9
    subplot(3,3,n)
    plot(mesh,V(:,n))
    title(['k = ',num2str(EV(n))])
    shading interp
    alpha(0.1)
    axis equal off
end

% Analytical solutions of eigenvalues for an arbitrary cube
ref = zeros(Neig^3,1);
l = 1;
for i = 1:Neig
    for j = 1:Neig
        for k = 1:Neig
            ref(l) = pi*sqrt( (i/L(1))^2 + (j/L(2))^2 +  (k/L(3))^2 );
            l = l+1;
        end
    end
end
ref = sort(ref);
ref = ref(1:Neig);

% Error
sol = EV(1:Neig);
[ref sol abs(sol-ref)./ref]


disp('Done, thanks for use.')
