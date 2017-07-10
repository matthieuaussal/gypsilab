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
%| Author(s)  : François Alouges - CMAP, Ecole polytechnique              |
%| Creation   : 06.07.17                                                  |
%| Last modif : 06.07.17                                                  |
%| Synopsis   : Test of the finite element method                         |
%+========================================================================+

% Cleaning
clear all
close all
clc

% Library path
addpath('../openFem')
addpath('../openHmx')
addpath('../openMsh')

errL2P1 = []; errH1P1 = errL2P1; errL2P2 = []; errH1P2 = errL2P2; h = [];
Uex = @(x) cos(x(:,1)).*cos(x(:,2))/3;
f = @(x) cos(x(:,1)).*cos(x(:,2));


for nbPts=10:3:80
    % Square mesh
    Omega = mshSquare(nbPts^2,[2*pi 2*pi]);
    Omega.gss = 7;
    h = [h, 2*pi/nbPts];
    %plot(Omega);
    
    
    %%% SOLVE LINEAR PROBLEM
    disp('~~~~~~~~~~~~~ SOLVE LINEAR PROBLEM ~~~~~~~~~~~~~')
    
    
    % Finite elements
    u = fem('P1'); v = fem('P1');
    u2 = fem('P2'); v2 = fem('P2');
    
    % Matrix and RHS
    K = integral(Omega,grad(u),grad(v)) + integral(Omega,u,v);
    K2 = integral(Omega,grad(u2),grad(v2)) + integral(Omega,u2,v2);
    F = integral(Omega, v, f);
    F2 = integral(Omega, v2, f);
    
    % Résolution
    uh = K\F;
    uh2 = K2\F2;

    % erreur en norme L2 et H1
    errL2P1 = [errL2P1, u.diff(Omega, uh, Uex, 'L2')];
    errH1P1 = [errH1P1, u.diff(Omega, uh, Uex,  'H1')];
    errL2P2 = [errL2P2, u2.diff(Omega, uh2, Uex,  'L2')];
    errH1P2 = [errH1P2, u2.diff(Omega, uh2, Uex,  'H1')];
end

% Plot the error graphs
figure(2)
subplot(1,2,1)
loglog(h,errL2P1,'b+-',h,errL2P2,'r+-',h,(h*10^-0.5).^2,'k--',h,(h*10^-0.5).^3,'k:')
legend('EF P1', 'EF P2','slope 2','slope 3')
title('error L2, CB Neumann')

xlabel('h')
subplot(1,2,2)
loglog(h,errH1P1,'b+-',h,errH1P2,'r+-',h,h,'k',h,(h*10^-0.5).^2,'k--')
legend('EF P1', 'EF P2','slope 1','slope 2')
title('error H1, CB Neumann')
xlabel('h')

disp('Done.')