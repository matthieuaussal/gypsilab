%|                                                                        |
%|                      THE GYPSILAB TOOLBOX - v0.20                      |
%|               www.cmap.polytechnique.fr/~aussal/gypsilab               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
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
%| Synopsis   : Triangles geometry and finite element validation          |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openMsh')
addpath('../openHmx')

% Parameters
Nvtx = 1e2;
tol  = 1e-3;


%%%%%%%%%%%%%%% GEOMETRY %%%%%%%%%%%%%%%
disp('=========== GEOMETRY ============')

% Create mesh disk
mesh = mshDisk(Nvtx,1);

% Graphical representation
plot(mesh); 
axis equal;
title('Mesh representation')
xlabel('X');   ylabel('Y');   zlabel('Z');
alpha(0.99)
hold on

% Centers
Xctr  = mesh.ctr;
meshc = msh(Xctr,(1:length(Xctr))');
if Nvtx < 1e3
    plot(meshc,'*b')
end

% Submeshing
Ip    = find(sum(Xctr>0,2)==2);
meshp = mesh.sub(Ip);
if Nvtx < 1e3
    plot(meshp,1)
    alpha(0.5)
end

% Edges
[edg2vtx,elt2edg] = mesh.edg;
tmp   = elt2edg(Ip,:);
tmp   = unique(tmp);
meshe = msh(mesh.vtx,edg2vtx(tmp,:));
if Nvtx < 1e3
    plot(meshe,'.-r')
end

% Free boundaries
bnd2vtx = mesh.bnd;
meshb   = msh(mesh.vtx,bnd2vtx);
if Nvtx < 1e3
    plot(meshb,'o-b')
end

% Surface
A = mesh.ndv;
abs(sum(A) - pi)/pi
P = meshb.ndv;
abs(sum(P) - 2*pi)/(2*pi)

% Normals
pts = mesh.ctr;
vct = mesh.nrm;
if Nvtx < 1e3
    quiver3(pts(:,1),pts(:,2),pts(:,3),vct(:,1),vct(:,2),vct(:,3),'k');
end

% Edges normals
Nu = mesh.nrmEdg;
if Nvtx < 1e3
    for i = 1:3
        pts = mesh.vtx(mesh.elt(Ip,i),:);
        vct = Nu{i}(Ip,:);
        quiver3(pts(:,1),pts(:,2),pts(:,3),vct(:,1),vct(:,2),vct(:,3),'g');
    end
end

% Surface quadrature
mesh.gss = 3;
[Xqud,W] = mesh.qud;
if Nvtx < 1e3
    plot3(Xqud(:,1),Xqud(:,2),Xqud(:,3),'xm')
    alpha(0.99)
end
abs(sum(W) - pi)/pi

% Circum quadrature
meshb.gss = 1;
[xqud,w] = meshb.qud;
if Nvtx < 1e3
    plot3(xqud(:,1),xqud(:,2),xqud(:,3),'xm')
    alpha(0.99)
end
abs(sum(w) - 2*pi)/(2*pi)

% Quadrature normals
pts = Xqud;
vct = mesh.nrmQud;
if Nvtx < 1e3
    quiver3(pts(:,1),pts(:,2),pts(:,3),vct(:,1),vct(:,2),vct(:,3),'m');
end

% Finite element dof
u    = fem('P1');
v    = fem('P1');
Xdof = v.dof(mesh);
if Nvtx < 1e3
    plot3(Xdof(:,1),Xdof(:,2),Xdof(:,3),'ob')
end


%%%%%%%%%%%%%%% SINGLE INTEGRATION %%%%%%%%%%%%%%%
disp('=========== SINGLE INTEGRATION ============')

% Numerical function
Id  = @(X) ones(size(X,1),1);

%-------------------------------------

% \int_{mesh(x)} f(x) dx 
I = integral(mesh,Id);
abs(I-pi)/pi

%-------------------------------------

% \int_{mesh(x)} f(x) psi(x) dx 
I = integral(mesh,Id,v);
abs(sum(I,2) - pi)/pi

% \int_{bnd(x)} f(x) psi(x) dx 
I = integral(meshb,Id,v);
P = v.sub(meshb,mesh);
I = I * P;
abs(sum(I,2) - 2*pi)/(2*pi)

%-------------------------------------

% \int_{mesh(x)} psi(x)' f(x)  dx 
I = integral(mesh,u,Id);
abs(sum(I,1) - pi)/pi

%-------------------------------------

% \int_{mesh(x)} psi(x)' psi(x) dx 
I = integral(mesh,u,v);
abs(sum(sum(I,1),2) - pi)/pi

% \int_{mesh(x)} grad(psi)(x)' grad(psi(x)) dx 
I = integral(mesh,grad(u),grad(v));
sum(sum(I,1),2)

%-------------------------------------

% \int_{mesh(x)} psi(x)' f(x) psi(x) dx 
I = integral(mesh,u,Id,v);
abs(sum(sum(I,1),2) - pi)/pi

% \int_{mesh(x)} grad(psi(x))' f(x) grad(psi(x)) dx 
I = integral(mesh,grad(u),Id,grad(v));
sum(sum(I,1),2)


%%%%%%%%%%%%%%% DOUBLE INTEGRATION %%%%%%%%%%%%%%%
disp('=========== DOUBLE INTEGRATION ============')

% Numerical function
Id  = @(X,Y) (X(:,1)==Y(:,1)) .* (X(:,2)==Y(:,2)) ;

% \int_{mesh(y)} f(x,y) psi(y) dy    
I = integral(mesh.qud,mesh,Id,v);
abs(sum(sum(I,1),2) - pi)/pi

%-------------------------------------

% \int_{mesh(x)} psi(x)' f(x,y) dx  
I = integral(mesh,mesh.qud,u,Id);
abs(sum(sum(I,1),2) - pi)/pi

%-------------------------------------

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dx 
I = integral(mesh,mesh,u,Id,v);

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' f(x,y) grad(psi(y)) dx 
I = integral(mesh,mesh,grad(u),Id,grad(v));
sum(sum(I,1),2)


%%%%%%%%%%%%%%% H-MATRIX INTEGRATION %%%%%%%%%%%%%%%
disp('=========== H-MATRIX INTEGRATION ============')

% Numerical function
G = @(X,Y) 1./(4*pi) * femGreenKernel(X,Y,'[exp(ikr)/r]',5);

% \int_{mesh(y)} G(x,y) psi(y) dy    
sol = integral(mesh.ctr,mesh,G,v,tol);
ref = integral(mesh.ctr,mesh,G,v);
norm(full(sol)-ref)./norm(ref)

%-------------------------------------

% \int_{mesh(x)} psi(x)' G(x,y) dx    
sol = integral(mesh,mesh.ctr,u,G,tol);
ref = integral(mesh,mesh.ctr,u,G);
norm(full(sol)-ref)./norm(ref)

%-------------------------------------

% \int_{mesh(x)} \int_{mesh(y)} psi(x)' G(x,y) psi(y) dx dy   
sol = integral(mesh,mesh,u,G,v,tol);
ref = integral(mesh,mesh,u,G,v);
norm(full(sol)-ref)./norm(ref)

% \int_{mesh(x)} \int_{mesh(y)} grad(psi(x))' G(x,y) grad(psi(y)) dx    
sol = integral(mesh,mesh,nxgrad(u),G,nxgrad(v),tol);
ref = integral(mesh,mesh,nxgrad(u),G,nxgrad(v));
norm(full(sol)-ref)./norm(ref)


disp('Done, thanks for use.')


