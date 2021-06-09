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
%|    #    |   FILE       : nrt2DiskUnbaffled.m                           |
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
name  = '2DiskUnbaffled'
N     = 200
tol   = 1e-3
typ   = 'P1'
gss   = 3
a1    = 0.05
a2    = 0.03
L     = 0.10
nk    = 64
kamax = 3*pi-0.2
v0    = 1

% Wave number discretization
K = linspace(1e-3,kamax,nk)/a1;

% Disk 1
disk = mshDisk(N,a1);
Xdwn = disk.vtx;
edge = disk.stp;
lgt  = edge(3);
Xup  = Xdwn;
Xup(:,3) = Xup(:,3) + lgt/10;
X    = [Xdwn;Xup];
tet = delaunayTriangulation(X(:,1),X(:,2),X(:,3));
[tri,vtx] = freeBoundary(tet);
disk1 = msh(vtx,tri);

% Disk 2
disk = mshDisk(N,a2);
Xdwn = disk.vtx;
edge = disk.stp;
Xup  = Xdwn;
Xup(:,3) = Xup(:,3) + edge(3);
X = [Xdwn;Xup];
X(:,3) = X(:,3) + L + lgt/10;  
tet = delaunayTriangulation(X(:,1),X(:,2),X(:,3));
[tri,vtx] = freeBoundary(tet);
disk2 = msh(vtx,tri);

% Check L
if (min(disk2.vtx(:,3)) - max(disk1.vtx(:,3)) ~= L)
    error('epaisseur non valide')
end

% Final mesh
mesh = union(disk1,disk2)

% Domain + finite element space
sigma = dom(mesh,gss);  % quadrature
u     = fem(mesh,typ);  % unknown
v     = fem(mesh,typ);  % test fct p1
w     = fem(mesh,'P0'); % test fct p0

% Quadrature on Sigma
[Xqud,Wqud] = sigma.qud;
Nqud        = sigma.qudNrm;
Mq2ddl      = u.uqm(sigma);    
Iqud1       = find(Xqud(:,3)<L/2);
Iqud2       = find(Xqud(:,3)>L/2);

% Mass matrices
Mp1  = integral(sigma,v,u);   % mass matrix P1P1
Mp10 = integral(sigma,v,w);   % mass matrix P1P0

% Speed
V0    = zeros(size(mesh.elt,1),1);
ctr   = mesh.ctr;
ctr1  = disk1.ctr;
z1min = min(ctr1(:,3));
z1max = max(ctr1(:,3));
V0(ctr(:,3)==z1min) = -v0;
V0(ctr(:,3)==z1max) = v0;

% Affichage
% figure(1)
% plot(mesh,V0)
% hold on
% plotNrm(mesh)
% plot(msh(Xqud(Iqud1,:)),'r');
% plot(msh(Xqud(Iqud2,:)),'y');
% quiver3(Xqud(:,1),Xqud(:,2),Xqud(:,3),Nqud(:,1),Nqud(:,2),Nqud(:,3),'k')
% axis equal
% % alpha(0.5)
% grid on
% colorbar 


%%% SOLVER %%%
disp('Solving...')

% Prepare impedance
Z11      = cell(size(K));
Z12      = cell(size(K));

% Parallel loop on k
tic
parfor l = 1:length(K)
    % Wave number
    k = K(l);
    
    % Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y|
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
    
    % Forces 
    % \int_sigma p(x)n(x) d_sx = sum_g w_g * p(x_g) * nz(x_g)
    Fz1 = sum(Wqud(Iqud1) .* Pqud(Iqud1) .* Nqud(Iqud1,3));
    Fz2 = sum(Wqud(Iqud2) .* Pqud(Iqud2) .* Nqud(Iqud2,3));
    
    % Impedances Zij = \frac{Fi}{si*vj}  (v1=v0, v2=0)
    Z11{l} = 1i*k * Fz1/(sum(disk1.ndv)*v0);
    Z12{l} = 1i*k * Fz2/(sum(disk1.ndv)*v0);
    
    % Affichage
%     figure
%     plot(mesh,real(pe))
%     axis equal
%     grid on
%     colorbar
end
toc

% Impedance finale avec conversion psi et surface
Z11 = cell2mat(Z11);
Z12 = cell2mat(Z12);

% Analytical impedance
load(['ref',name,'.mat']);


%%% COMPARAISON
% Real part
figure(2)
subplot(2,1,1)
plot(Zref.ka,real(Zref.z),'-k')
hold on
plot(K*a1,real(Z11),'xk')
grid on 
legend({'Analytic','Numeric'},'Location','NorthWest')
xlabel('Non-dimensional frequency ka')
ylabel('Re(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Imaginary part
subplot(2,1,2)
plot(Zref.ka,imag(Zref.z),'-k')
hold on
plot(K*a1,-imag(Z11),'xk')
grid on 
legend({'Analytic','Numeric'},'Location','SouthWest')
xlabel('Non-dimensional frequency ka')
ylabel('Im(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Export figure
figure(2)
print('-dpng','-r300',['fig',name])



disp('~~> Michto gypsilab !')




