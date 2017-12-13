function [X,elt2dof,col] = femDof(fe)
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
%|    #    |   FILE       : femDof.m                                      |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Define degrees of freedom                     |
%|  `---'  |                                                              |
%+========================================================================+

% Initialization
mesh    = fe.msh;
c       = unique(mesh.col);
Ic      = cell(length(c),1);
X       = cell(length(c),1);
elt2dof = cell(length(c),1);
col     = cell(length(c),1);
Nx      = 0;

% Loop on colors
for i = 1:length(c)
    % Submeshing with uniform colour
    Ic{i}  = find(mesh.col==c(i));
    fe.msh = mesh.sub(Ic{i});
    
    % Lagrange order 0, constant by element
    if strcmp(fe.typ,'P0')
        X{i}       = fe.msh.ctr;
        elt2dof{i} = (1:size(fe.msh.elt,1))';
        
    % Lagrange order 1, piecewise linear by element
    elseif strcmp(fe.typ,'P1')
        X{i}       = fe.msh.vtx;
        elt2dof{i} = fe.msh.elt;
        
    % Lagrange order 2, piecewise quadratic by element
    elseif strcmp(fe.typ,'P2')
        [edge,elt2edg] = fe.msh.edg;
        X{i}           = [fe.msh.vtx ; edge.ctr];
        elt2dof{i}     = [fe.msh.elt , elt2edg + size(fe.msh.vtx,1)];
        
    % NED linear elements
    elseif strcmp(fe.typ,'NED')
        [edge,elt2dof{i}] = fe.msh.edg;
        X{i}              = edge.ctr;
        
    % RWG linear elements for surfacic mesh
    elseif strcmp(fe.typ,'RWG')
        % Triangular mesh
        if size(fe.msh.elt,2)==3
            [edge,elt2dof{i}] = fe.msh.edg;
            X{i}              = edge.ctr;
        
        % Tetrahedral mesh
        elseif size(fe.msh.elt,2)==4
            [face,elt2dof{i}] = fe.msh.fce;
            X{i}              = face.ctr;
        end
        
    % Others
    else
        error('femDof.m : unavailable case')
    end
    
    % Colors
    col{i} = c(i) * ones(size(X{i},1),1);
    
    % Incrementation
    elt2dof{i} = Nx + elt2dof{i};
    Nx         = Nx + size(X{i},1);
end

% Vectorial format
X       = cell2mat(X);
elt2dof = cell2mat(elt2dof);
col     = cell2mat(col);
Ic      = cell2mat(Ic);

% Security for lost elements
if (norm(sort(Ic) - (1:size(elt2dof,1))','inf') > 1e-12)
   error('femDom.m : unavailable case') 
end

% Reordering elements
elt2dof(Ic,:) = elt2dof;
end

