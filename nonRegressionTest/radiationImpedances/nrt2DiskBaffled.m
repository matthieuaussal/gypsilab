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
%|    #    |   FILE       : nrt2DiskBaffled.m                             |
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
name  = '2DiskBaffled'
N     = 200
gss   = 3
a     = 0.05
kamax = 3*pi-0.2
v0    = 1
nk    = 64
L     = 0.1
nimag = 100 

% Wave number discretization
K = linspace(1e-3,kamax,nk)/a;

% Prepare mesh, domain and finite element space
mesh  = mshDisk(N,a)    % mesh
sigma = dom(mesh,gss);  % quadrature 
u     = fem(mesh,'P1'); % unknown P1

% Quadrature on Sigma
sigmaDisk   = dom(mesh,gss);  % quadrature 
[Xqud,Wqud] = sigmaDisk.qud;
area        = sum(mesh.ndv);

% Regularisation 
Greg = full(sum(regularize(Xqud,sigma,'[1/r]',u),2));

% Affichage
figure(1)
set(gcf,'Color','w')
plot(mesh,'w')
hold on
plot3(Xqud(:,1),Xqud(:,2),Xqud(:,3),'.k')
axis equal off
grid on


%%% SOLVER %%%
disp('Solving...')

% Prepare impedance
Z = cell(nimag+1,1);

% Loop on kernel
tic
parfor nc = 1:length(Z)
    % Image
    i = nc-1;
    
    % Distances for Gauss quadrature
    ind   = 1:length(Wqud);
    [I,J] = ndgrid(ind,ind);
    Rxy   = sqrt( ...
        (Xqud(I,1)-Xqud(J,1)).^2 + ...
        (Xqud(I,2)-Xqud(J,2)).^2 + ...
        (i*2*L).^2 );
    Rxy   = reshape(Rxy,length(Wqud),length(Wqud));
    
    % Loop on wave number
    tmp = zeros(size(K));
    for j = 1:length(K)
        % Wave number
        k = K(j);
        
        % Green kernel
        Gxy = exp(1i*k*Rxy)./Rxy;
        
        % Pressure
        if (i==0)
            Gxy(Rxy<1e-12) = 0 + 1i*k;
            Pqud = -2/(4*pi) * (Gxy*Wqud + Greg) * v0;
        else
            Pqud = -2/(4*pi) * (Gxy*Wqud) * v0;
        end
        
        % Force
        % \int_sigma p(x) n(x) d_sx = sum_g w_g * p(x_g)
        Fz = Wqud' * Pqud;
        
        % Impedance Z = \frac{Fz}{s*v}
        if (i==0)
            tmp(j) = 1i*k * Fz/(area*v0);
        else % for symetrie
            tmp(j) = 2i*k * Fz/(area*v0);
        end
    end
    
    % Output
    Z{nc} = tmp;
end
toc     

% Numerical impedance
Zsol = sum(cell2mat(Z),1);

% Analytical impedance
load(['ref',name,'.mat']);
I = find(Zref.ka<kamax);
Zref.ka = Zref.ka(I);
Zref.z  = Zref.z(I);

% Tube impedance
Ztub = -1./(tan(L/a.*Zref.ka));
Ztub = min(max(Ztub,-10),10);



%%% COMPARAISON
% Real part
figure(2)
subplot(2,1,1)
plot(Zref.ka,real(Zref.z),'-k')
hold on
plot(K*a,real(Zsol),'xk')
grid on 
legend({'Analytic','Numeric'})
xlabel('Non-dimensional frequency ka')
ylabel('Re(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Imaginary part
subplot(2,1,2)
plot(Zref.ka,imag(Zref.z),'-k')
hold on
plot(Zref.ka,Ztub,':k')
plot(K*a,-imag(Zsol),'xk')
grid on 
legend({'Analytic','Tube','Numeric'})
xlabel('Non-dimensional frequency ka')
ylabel('Im(Z/Z_0)')
set(gca,'XTick',0:pi/2:3*pi)
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi','5\pi/2','3\pi'})

% Export figure
figure(1)
print('-dpng','-r300',['fig','DiskMesh'])
figure(2)
print('-dpng','-r300',['fig',name])



disp('~~> Michto gypsilab !')


