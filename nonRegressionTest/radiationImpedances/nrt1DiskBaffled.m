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
%|    #    |   FILE       : nrt1DiskBaffled.m                             |
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
name  = '1DiskBaffled'
N     = 200
gss   = 3
a     = 0.05
kamax = 3*pi-0.2
v0    = 1
nk    = 64

% Wave number discretization
K = linspace(1e-3,kamax,nk)/a;

% Prepare 2D mesh, domain and finite element space
mesh  = mshDisk(N,a)    % mesh
sigma = dom(mesh,gss);  % quadrature 
u     = fem(mesh,'P1'); % unknown P1
v     = fem(mesh,'P1'); % test P1
w     = fem(mesh,'P0'); % test P0

% Quadrature on Sigma
[Xqud,Wqud] = sigma.qud;
Mq2ddl      = u.uqm(sigma);    

% Mass matrices
Mp1  = integral(sigma,v,u);   % mass matrix P1P1
Mp10 = integral(sigma,v,w);   % mass matrix P1P0

% Speed
V0    = zeros(size(mesh.elt,1),1);
V0(:) = v0;

% Visu
figure(1)
plot(mesh)
hold on
plotNrm(mesh)
plot(msh(sigma.qud),'r')
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
    Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);
    
    % Single layer operators 
    Lk_SS = 1/(4*pi) .* integral(sigma,sigma,v,Gxy,u);
    Lk_SS = Lk_SS + 1/(4*pi) .* regularize(sigma,sigma,v,'[1/r]',u);
    
    % Normal derivative
    v_S = Mp1\(Mp10*V0);
    
    % Rayleigh integral on surface
    phi_s = -2 * (Mp1\(Lk_SS * v_S));
    
    % Disk force \int_sigma p(x) n(x) d_sx = sum_g w_g * p(x_g)
    Pqud = Mq2ddl * phi_s;
    Fz   = sum(Wqud .* Pqud);
    
    % Impedance Z = \frac{Fz}{s*v}
    Z{l} = 1i*k * Fz/(sum(mesh.ndv)*v0);
end
toc

% Numerical impedance
Zsol = cell2mat(Z);

% Analytical impedance
load(['ref',name,'.mat']);


%%% COMPARAISON
% Real part
figure(2)
subplot(2,1,1)
plot(Zref.ka,real(Zref.z),'-k')
hold on
plot(K*a,real(Zsol),'xk')
grid on 
legend({'Analytic','Numeric'},'Location','SouthEast')
xlabel('Non-dimensional frequency ka')
ylabel('Re(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Imaginary part
subplot(2,1,2)
plot(Zref.ka,imag(Zref.z),'-k')
hold on
plot(K*a,-imag(Zsol),'xk')
grid on 
legend({'Analytic','Numeric'})
xlabel('Non-dimensional frequency ka')
ylabel('Im(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Export figure
figure(2)
print('-dpng','-r300',['fig',name])



disp('~~> Michto gypsilab !')




