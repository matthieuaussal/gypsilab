function err = femDiff(fe, mesh, Uh, Uex, type )
%+========================================================================+
%|                                                                        |
%|           OPENMSH, MESH MANAGEMENT AND NUMERICAL QUADRATURE            |
%|              openMsh is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
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
%| Author(s)  : François Alouges - CMAP, Ecole polytechnique              |
%| Creation   : 06.07.17                                                  |
%| Last modif : 06.07.17                                                  |
%| Synopsis   : L2 and H1 errors                                          |
%+========================================================================+

% Mesh integration data
[X,Wx] = mesh.qud;

% Finite element matrix
Mu = fe.dqm(mesh);

% Function to be applied
Uexact = Uex(X);
if (size(Uexact,2) > 1)
    error('mshIntegral.m : unavailable case')
end
    
Uapp = Mu * Uh;
switch type
    case 'L2'
        err = sqrt(sum(Wx.*(Uexact-Uapp).^2));
    case 'H1'
        eps = 1e-6;
        gradef = grad(fe);
        Gu = gradef.dqm(mesh);
        DxUapp = Gu{1} * Uh;
        DyUapp = Gu{2} * Uh;
        DzUapp = Gu{3} * Uh;
        n1 = size(X,1);
        eps1 = eps*ones(n1,1)*[1,0,0];
        eps2 = eps*ones(n1,1)*[0,1,0];
        eps3 = eps*ones(n1,1)*[0,0,1];
        DxUexact = (Uex(X+eps1)-Uex(X-eps1))/(2*eps);
        DyUexact = (Uex(X+eps2)-Uex(X-eps2))/(2*eps);
        DzUexact = (Uex(X+eps3)-Uex(X-eps3))/(2*eps);
        err = sqrt(sum(Wx.*((Uexact-Uapp).^2+(DxUexact-DxUapp).^2+(DyUexact-DyUapp).^2+(DzUexact-DzUapp).^2)));
    otherwise
        error('Unknown error type. Known types are ''L2'' or ''H1''.');
end
end

