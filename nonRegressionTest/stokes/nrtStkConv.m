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
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : nrtStokesGTradGal.m                           |
%|    #    |   VERSION    : 0.42                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 31.10.2018                                    |
%|  / 0 \  |   LAST MODIF : 31.10.2018                                    |
%| ( === ) |   SYNOPSIS   : Stoke Galerkin of an ovoid using BEM with     |
%|  `---'  |                stokeslet G and stresslet T                   |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../../openDom')
addpath('../../openFem')
addpath('../../openMsh')

% Parameters
N   = 2.^(5:7);
x0  = [1 2 0.5];
err = zeros(size(N));
h   = zeros(size(N));

% Loop over mesh
for n = 1:length(N)
    % Mesh unit sphere
    mesh = mshSphere(N(n),1);
    
    % Sphere to ellipsoid
    mesh.vtx(:,1) = 5 * mesh.vtx(:,1);
    mesh.vtx(:,2) = 3 * mesh.vtx(:,2);
    mesh.vtx(:,3) = 2 * mesh.vtx(:,3);
    
    % Quadrature
    gamma = dom(mesh,3);
    
    % Finite element
    phi = fem(mesh,'P1');
    unk = phi.unk;
    un  = ones(length(phi),1);

    % Mass matrix
    M = integral(gamma,phi,phi);
    
    % Graphical rep
%     figure(n)
%     plot(mesh)
%     hold on
%     plotNrm(mesh)
%     plot(gamma)
%     plot(phi,'og')
%     plot3(x0(1),x0(1),x0(1),'*r')
%     axis equal
%     alpha(0.5)
    
    % Initalization
    G      = cell(3,3);
    T      = cell(3,3);
    C      = cell(3,3);
    I      = cell(3,3);
    mu     = cell(3,1);
    lambda = cell(3,1);
    
    % Loop for each coordinate
    for i = 1:3
        for j = 1:3
            % Using r = x - y
            % Single layer : G = \int_gamma \int_gamma  1/(8pi) (\delta_ij/r + r_i*r_j/|r|^3)
            name   = ['[ij/r+rirj/r^3]',num2str(i),num2str(j)];
            green  = @(X,Y) 1/(8*pi) .* femGreenKernel(X,Y,name,[]);
            G{i,j} = integral(gamma,gamma,phi,green,phi);
            
            % Regularization
            G{i,j} = G{i,j} + 1/(8*pi) .* regularize(gamma,gamma,phi,name,phi);
            
            % Double layer  : T = \int_gamma \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
            % Double lumped : P = \int_gamma -6/(8pi) (r_i*r_j*(r.n)/|r|^5)
            T{i,j} = 0;
            P      = 0;
            for k = 1:3
                name   = ['[rirjrk/r^5]',num2str(i),num2str(j),num2str(k)];
                green  = @(X,Y) -6/(8*pi) .* femGreenKernel(X,Y,name,[]);
                T{i,j} = T{i,j} + integral(gamma,gamma,phi,green,ntimes(phi,k));
                P      = P + integral(gamma.qud,gamma,green,ntimes(phi,k)) * un;
            end
            
            % Correction
            C{i,j} = integral(gamma,phi,@(X)P,phi);
            
            % Identity matrix
            if (i == j)
                I{i,j} = M;
            else
                I{i,j} = sparse(size(M,1),size(M,2));
            end
        end
        
        % mu = [u] = u_int - u_ext = - Gi1 = - Gi1(x0,y)
        name  = ['[ij/r+rirj/r^3]',num2str(i),num2str(1)];
        green = @(Y) 1/(8*pi) .* femGreenKernel(x0,Y,name,[]);
        mu{i} = - integral(gamma,phi,green);
        
        % lambda = [sigma] = - Ti1
        lambda{i} = 0;
        for k = 1:3
            name      = ['[rirjrk/r^5]',num2str(i),num2str(1),num2str(k)];
            green     = @(Y) -6/(8*pi) .* femGreenKernel(x0,Y,name,[]);
            lambda{i} = lambda{i} - integral(gamma,ntimes(phi,k),green);
        end
    end
    
    % Convert cell to full matrix
    G      = cell2mat(G);
    T      = cell2mat(T);
    C      = cell2mat(C);
    I      = cell2mat(I);
    mu     = cell2mat(mu);
    lambda = cell2mat(lambda);
    
    % Analytic solution : Gi1
    ref = - mu;
        
    % Stokes radiation with correction
    sol = -G*(I\lambda) + T*(I\mu) - C*(I\mu);
%     sol = -G*(I\lambda) + T*(I\mu) - 0.5*mu;

    % Relative error L2 and inf
    tmp    = ref-sol;
    err(n) = sqrt(tmp'*(I\tmp));
    h(n)   = mean(mesh.stp);
    norm(tmp)./norm(ref)
end

% Graphical representation
figure(1)
hold on
loglog(h,err,'-ok')
grid on
% legend({'G & T','G','Id/r','0'})


disp('~~> Michto gypsilab !')




