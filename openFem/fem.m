classdef fem
%+========================================================================+
%|                                                                        |
%|                  OPENFEM, FINITE AND BOUNDARY ELEMENT                  |
%|              openFem is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openFem,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Finite element class definition                           |
%+========================================================================+

properties
    typ = [];         % FINITE ELEMENT TYPE (P0, P1)
    opr = [];         % OPERATOR APPLED TO FINITE ELEMENT
end

methods
    % CONSTRUCTOR
    function fe = fem(str)
        fe.typ = str;
        fe.opr = '[psi]';
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
    
    % NORMAL DOT BASIS FUNCTION
    function fe = ndot(varargin)
        fe = varargin{1};
        if (nargin == 1)
            fe.opr = 'n*[psi]';
        else
            fe.opr = ['n*[psi]',num2str(varargin{2})];
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
    
    % DEGREES OF FREEDOM
    function [X,elt2dof] = dof(fe,mesh)
        [X,elt2dof] = femDof(fe,mesh);
    end
    
    % L2 AND H1 ERRORS
    function err = diff(fe, mesh, sol, ref, type )
        err = femDiff(fe, mesh, sol, ref, type );
    end
    
    % DOF TO QUADRATURE MATRIX
    function M = dqm(fe,mesh)
        if strcmp(fe.typ(1),'P')
            M = femLagrangePn(fe,mesh);
        else
            erro('fem.m : unavailable case')
        end
    end      
    
    % NUMERICAL REGULARIZAION
    function S = regularize(varargin)
        S = femRegularize(varargin);
    end
            
    % SUBDOMAIN TRANSFERT MATRIX
    function M = sub(u,mesh1,mesh2)
        dof1      = u.dof(mesh1);
        dof2      = u.dof(mesh2);
        [~,I1,I2] = intersect(dof1,dof2,'rows');
        M         = sparse(I1,I2,1,size(dof1,1),size(dof2,1));
    end
    
    % DIRICHLET ELIMINATION MATRIX
    function M = dir(u,mesh1,mesh2)
        dof1      = u.dof(mesh1);
        dof2      = u.dof(mesh2);
        [~,nodir] = setdiff(dof1,dof2,'rows');
        M         = sparse(nodir,1:length(nodir),1,size(dof1,1),length(nodir));
    end
end
end
