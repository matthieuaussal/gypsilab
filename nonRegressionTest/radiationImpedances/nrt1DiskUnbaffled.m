%+========================================================================+
%|                                                                        |
%|            This script uses the GYPSILAB toolbox for Matlab            |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrt1DiskUnbaffled.m                           |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & Joel Bensoam                |
%|  ( # )  |   CREATION   : 04/06/2021                                    |
%|  / 0 \  |   LAST MODIF :                                               |
%| ( === ) |   SYNOPSIS   :                                               |
%|  `---'  |                                                              |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
run('../../addpathGypsilab.m')

% Parametres
name  = '1DiskUnbaffled'
N     = 200
tol   = 1e-3
typ   = 'P1'
gss   = 3
a     = 0.05
nk    = 64
kamax = 3*pi-0.2
v0    = 1

% Wave number discretization
K = linspace(1e-3,kamax,nk)/a;

% Disk mesh (3D)
disk = mshDisk(N,a);
Xdwn = disk.vtx;
edge = disk.stp;
lgt  = edge(3);
Xup  = Xdwn;
Xup(:,3) = Xup(:,3) + lgt/10;
X    = [Xdwn;Xup];
tet = delaunayTriangulation(X(:,1),X(:,2),X(:,3));
[tri,vtx] = freeBoundary(tet);
mesh = msh(vtx,tri)

% Domain + finite element space
sigma = dom(mesh,gss);  % quadrature
u     = fem(mesh,'P1'); % unknown P1
v     = fem(mesh,'P1'); % test P1
w     = fem(mesh,'P0'); % test P0

% Quadrature on Sigma
[Xqud,Wqud] = sigma.qud;
Nqud        = sigma.qudNrm;
Mq2ddl      = u.uqm(sigma);    

% Mass matrices
Mp1  = integral(sigma,v,u);   % mass matrix P1P1
Mp10 = integral(sigma,v,w);   % mass matrix P1P0

% Speed
V0    = zeros(size(mesh.elt,1),1);
ctr   = mesh.ctr;
z1min = min(ctr(:,3));
z1max = max(ctr(:,3));
V0(ctr(:,3)==z1min) = -v0;
V0(ctr(:,3)==z1max) = v0;

% Visu
figure(1)
plot(mesh,V0)
hold on
plotNrm(mesh)
plot(msh(Xqud),'r');
quiver3(Xqud(:,1),Xqud(:,2),Xqud(:,3),Nqud(:,1),Nqud(:,2),Nqud(:,3),'k')
axis equal
grid on
colorbar 


%%% SOLVER %%%
disp('Solving...')

% Prepare impedance
Z = cell(size(K));

% Parallel loop on k
tic
parfor l = 1:length(K)
    % Wave number
    k = K(l);
    
    % Green kernel function G(x,y) = exp(ik|x-y|)/|x-y|
    Gxy      = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
    gradyGxy = {
        @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k) , ...
        @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k) , ...
        @(X,Y) femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k) };
    
    % Double layer operator
    % D = 1/(4pi) \int_\sigma_x \int_sigma_y v(x) dn_y(G(x,y)) u(y) ds_x ds_y
    D  = 1/(4*pi) .* integral(sigma,sigma,v,gradyGxy,ntimes(u),tol);
    Dr = 1/(4*pi) .* regularize(sigma,sigma,v,'grady[1/r]',ntimes(u));
    D  = D + Dr;
    
    % Hypersingular operator
    % H = 1/(4pi) \int_\sigma_x \int_sigma_y v(x) dn__x(dn_y(G(x,y))) u(y) ds_x ds_y
    % H = 1/(4pi) * (k^2 * \int_\sigma_x \int_sigma_y n.v(x) G(x,y) n.u(y) ds_x ds_y
    % - \int_\sigma_x \int_sigma_y nxgrad(v(x)) G(x,y) nxgrad(u(y)) ds_x ds_y)
    N  = 1/(4*pi) .* (k^2 * integral(sigma,sigma,ntimes(v),Gxy,ntimes(u),tol) ...
        - integral(sigma,sigma,nxgrad(v),Gxy,nxgrad(u),tol));
    Nr = 1/(4*pi) .* (k^2 * regularize(sigma,sigma,ntimes(v),'[1/r]',ntimes(u)) ...
        - regularize(sigma,sigma,nxgrad(v),'[1/r]',nxgrad(u)));
    N  = N + Nr;
    
    % Solve linear system : [-N]mu = dnPe
    mu = (-N)\(Mp10*V0);
    
    % Boundary pressure : Mp1 p = [-Id/2 - D]mu
    pe = Mp1\((-Mp1/2-D)*mu);

    % \int_sigma p(x) d_sx = sum_g w_g * p(x_g)
    Pqud = Mq2ddl * pe;
    
    % Disk force \int_sigma p(x)n(x) d_sx = sum_g w_g * p(x_g) * nz(x_g)
    Fz1 = sum(Wqud .* Pqud .* Nqud(:,3));
    
    % Impedance Zij = \frac{Fi}{si*vj}  (v1=v0, v2=0)
    Z{l} = 1i*k * Fz1/(sum(mesh.ndv)*v0);
end
toc

% Numerical impedance
Zsol = cell2mat(Z);

% Analytical impedance
load(['ref',name,'.mat']);

% Mellow impedance
load(['ref',name,'Mellow.mat']);


%%% COMPARAISON
% Real part
figure(2)
subplot(2,1,1)
plot(Zref.ka,real(Zref.z),'--k')
hold on
plot(Zmel.ka,real(Zmel.z),':k')
plot(K*a,real(Zsol),'xk')
grid on 
legend({'Analytic','Mellow','Numeric'},'Location','SouthEast')
xlabel('Non-dimensional frequency ka')
ylabel('Re(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Imaginary part
subplot(2,1,2)
plot(Zref.ka,imag(Zref.z),'--k')
hold on
plot(Zmel.ka,imag(Zmel.z),':k')
plot(K*a,-imag(Zsol),'xk')
grid on 
legend({'Analytic','Mellow','Numeric'})
xlabel('Non-dimensional frequency ka')
ylabel('Im(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Export figure
figure(2)
print('-dpng','-r300',['fig',name])



disp('~~> Michto gypsilab !')



