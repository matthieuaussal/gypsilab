classdef fem
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
%|    #    |   FILE       : fem.m                                         |
%|    #    |   VERSION    : 0.32                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Finite element class definition               |
%|  `---'  |                                                              |
%+========================================================================+

properties
    typ = [];         % FINITE ELEMENT TYPE (P0, P1, P2, RWG)
    opr = [];         % OPERATOR APPLIED TO FINITE ELEMENT
    msh = [];         % FINITE ELEMENT SPACE 
    dir = [];         % MESH FOR DIRICHLET CONDITION 
    ctn = [];         % COLOURS FOR CONTINUITY
    jct = [];         % COLOURS FOR JUNCTIONS
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function fe = fem(mesh,str)
        fe.typ = str;
        fe.opr = '[psi]';
        fe.msh = mesh;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plot(varargin)
        fe  = varargin{1};
        spc = 'ob';
        if (nargin == 2)
            spc = varargin{2};
        end
        X = fe.dof;
        plot3(X(:,1),X(:,2),X(:,3),spc)
    end
    
    function plotUnk(varargin)
        fe  = varargin{1};
        spc = 'ob';
        if (nargin == 2)
            spc = varargin{2};
        end
        X = fe.unk;
        plot3(X(:,1),X(:,2),X(:,3),spc)
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LENGTH
    function s = length(fe)
        s = size(fe.dof,1);
    end
    
    % SIZE
    function s = size(fe)
        s = size(fe.dof);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CURL OF THE BASIS FUNCTION
    function fe = curl(fe)
        fe.opr = 'curl[psi]';
    end
    
    % DIVERGENCE OF THE BASIS FUNCTION
    function fe = div(fe)
        fe.opr = 'div[psi]';
    end
    
    % GRADIENT OF THE BASIS FUNCTION
    function fe = grad(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'grad[psi]';
        else
            fe.opr = ['grad[psi]',num2str(varargin{2})];
        end
    end    
    
    % NORMAL TIMES BASIS FUNCTION
    function fe = ntimes(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'n*[psi]';
        else
            fe.opr = ['n*[psi]',num2str(varargin{2})];
        end
    end
    
    % NORMAL WEDGE BASIS FUNCTION
    function fe = nx(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'nx[psi]';
        else
            fe.opr = ['nx[psi]',num2str(varargin{2})];
        end
    end
        
    % NORMAL CROSS GRADIENT OF THE BASIS FUNCTION
    function fe = nxgrad(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'nxgrad[psi]';
        else
            fe.opr = ['nxgrad[psi]',num2str(varargin{2})];
        end
    end
    
    % DIVERGENCE NORMAL CROSS OF THE BASIS FUNCTION
    function fe = divnx(fe)
        fe.opr = 'curl[psi]';
    end    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEGREES OF FREEDOM %%%%%%%%%%%%%%%%%%%%%%
    % DEGREES OF FREEDOM
    function [X,elt2dof,col] = dof(fe)
        [X,elt2dof,col] = femDof(fe);
    end
        
    % DOF TO QUADRATURE MATRIX
    function M = dqm(fe,domain)
        [~,P] = fe.unk;    
        if (size(fe.msh.elt,2) == size(domain.msh.elt,2) + 1)
            bound  = fe.msh.bnd;
            P      = restriction(fe,bound) * P;
            fe.msh = bound;
        end        
        if strcmp(fe.typ(1),'P')
            M = femLagrangePn(fe,domain);
        elseif strcmp(fe.typ,'RWG')
            M = femRaoWiltonGlisson(fe,domain);
        elseif strcmp(fe.typ,'NED')
            M = femNedelec(fe,domain);
        else
            error('fem.m : unavailable case')
        end
        if iscell(M)
            M{1} = M{1} * P;
            M{2} = M{2} * P;
            M{3} = M{3} * P;
        else
            M = M * P;
        end
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNKNOWNS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIRICHLET
    function fe = dirichlet(fe,mesh)
       fe.dir = mesh; 
    end
    
    % CONTINUITY
    function fe = continuity(fe,col)
       fe.ctn = [fe.ctn ; col ]; 
    end
    
    % JUNCTION
    function fe = junction(fe,c,a)
       fe.jct = [c a]; 
    end
    
    % UNKNOWNS AND REDUCTION MATRIX
    function [X,P] = unk(fe)
        [X,P] = femUnknown(fe);
    end
    
    % UNKNOWNS DATA TO VERTEX
    function I = feval(v,f,mesh)
        mesh.col(:) = 0;
        if (size(mesh.elt,2) == 4)
            gss = 4;
        else
            gss = 3;
        end        
        domain = dom(mesh,gss);
        u      = fem(mesh,'P1');
        M      = integral(domain,u,u);
        Fv     = integral(domain,u,v);
        if iscell(Fv)
            I{1} = M \ (Fv{1} * f);
            I{2} = M \ (Fv{2} * f);
            I{3} = M \ (Fv{3} * f);
        else
            I = M \ (Fv * f);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% USEFULL MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%
    % RESTRICTION MATRIX
    function M = restriction(u,mesh)
        v         = fem(mesh,u.typ);
        [~,I1,I2] = intersect(v.dof,u.dof,'rows');
        M         = sparse(I1,I2,1,size(v.dof,1),size(u.dof,1));
    end

    % ELIMINATION MATRIX
    function M = elimination(u,mesh)
        v         = fem(mesh,u.typ);
        [~,nodir] = setdiff(u.dof,v.dof,'rows');
        M         = sparse(nodir,1:length(nodir),1,size(u.dof,1),length(nodir));
    end
end
end
