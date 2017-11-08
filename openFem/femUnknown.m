function [X,P] = femUnknown(fe)
%+========================================================================+
%|                                                                        |
%|              OPENFEM - LIBRARY FOR FINITE ELEMENT METHOD               |
%|           openFem is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal & Francois Alouges (c) 2015-2017.          |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%|             francois.alouges@polytechnique.edu                         |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : femUnknown.m                                  |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2017                                    |
%| ( === ) |   SYNOPSIS   : Unknowns and reduction matrix for constrained |
%|  `---'  |                finite elements                               |
%+========================================================================+

% Initialization with dofs
[X,~,C] = dof(fe);
P       = speye(size(X,1));

% Dirichlet
if ~isempty(fe.dir)
    fed = fem(fe.dir,fe.typ);
    I   = find(~ismember(X,fed.dof,'rows'));
    M   = sparse(I,1:length(I),1,size(X,1),length(I));
    P   = P * M;
    X   = X(I,:);
    C   = C(I,:);
end

% Continuity
if ~isempty(fe.ctn)
    for i = 1:size(fe.ctn,1)
        R                   = rand(size(X));
        R(C==fe.ctn(i,1),:) = 0;
        R(C==fe.ctn(i,2),:) = 0;
        [~,I,J] = unique(X+R,'rows','stable');
        M = sparse(1:length(J),J,1,length(J),max(J));
        P = P * M;
        X = X(I,:);
        C = C(I,:);
    end
end

% Junctions
if ~isempty(fe.jct)
    % Indices for each colors
    I1 = find(C==fe.jct(1));
    I2 = find(C==fe.jct(2));
    I3 = find(C==fe.jct(3));
    
    % Indices for intersection between each colors
    [~,I31] = intersect(X(I3,:),X(I1,:),'rows');
    [~,I32] = intersect(X(I3,:),X(I2,:),'rows');
    I3      = I3(intersect(I31,I32));
    [~,I13] = intersect(X(I1,:),X(I3,:),'rows');
    I1      = I1(I13);
    [~,I23] = intersect(X(I2,:),X(I3,:),'rows');
    I2      = I2(I23);
    
    % Final unknowns
    I = setdiff((1:size(X,1))',I3);
    
    % Reduction matrix 
    idx = [I;I3;I3];
    jdx = [I;I1;I2];
    val = [ones(size(I)) ; ...
        - fe.jct(4)./fe.jct(6) * ones(size(I3));...
        - fe.jct(5)./fe.jct(6) * ones(size(I3))];
    M   = sparse(idx,jdx,val,size(X,1),size(X,1));
    M   = M(:,I);
    
    % Update
    P = P * M;
    X = X(I,:);
    C = C(I,:);
%     figure
%     spy(M==1,'b')
%     hold on
%     spy(M==-1,'r')
end
end
