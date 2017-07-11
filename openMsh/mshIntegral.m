function I = mshIntegral(data)
%+========================================================================+
%|                                                                        |
%|           OPENMSH, MESH MANAGEMENT AND NUMERICAL QUADRATURE            |
%|              openMsh is part of GYPSYLAB toolbox - v0.20               |
%|                                                                        |
%| Copyright (c) 20015-2017, Ecole polytechnique, all rights reserved.    |
%| Licence Creative Commons BY-NC-SA 4.0, Attribution, NonCommercial and  |
%| ShareAlike (see http://creativecommons.org/licenses/by-nc-sa/4.0/).    |
%| This software is the property from Centre de Mathematiques Appliquees  |
%| de l'Ecole polytechnique, route de Saclay, 91128 Palaiseau, France.    |
%|                                                            _   _   _   |
%| Please acknowledge the GYPSILAB toolbox in programs       | | | | | |  |
%| or publications in which you use the code. For openMsh,    \ \| |/ /   |
%| we suggest as reference :                                   \ | | /    |
%| [1] : www.cmap.polytechnique.fr/~aussal/gypsilab             \   /     |
%|                                                               | |      |
%|_______________________________________________________________|_|______|
%| Author(s)  : Matthieu Aussal - CMAP, Ecole polytechnique               |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Numerical integation                                      |
%+========================================================================+


%%% NUMERICAL INTEGRATION --> \int_{mesh(x)} f(x) dx 
if (length(data) == 2)
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    
    % Function to be applied
    F  = data{2};
    Fx = F(X);
    if (size(Fx,2) > 1)
        error('mshIntegral.m : unavailable case')
    end
    
    % Integration
    I = Wx' * Fx;
  
    
%%% FINITE ELEMENT INTEGRATION --> \int_{mesh(x)} f(x) psi(x) dx 
elseif (length(data) == 3) && isa(data{2},'function_handle') && isa(data{3},'fem')
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Wx     = spdiags(Wx,0,length(Wx),length(Wx));
    
    % Function to be applied
    F  = data{2};
    Fx = F(X);
    if (size(Fx,2) > 1)
        error('mshIntegral.m : unavailable case')
    end
    
    % Finite element matrix
    v  = data{3};
    Mv = v.dqm(Xmsh);
    
    % Integration
    I = Fx.' * Wx * Mv;
    

%%% FINITE ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x)  dx 
elseif (length(data) == 3) && isa(data{2},'fem') && isa(data{3},'function_handle')
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Wx     = spdiags(Wx,0,length(Wx),length(Wx));
    
    % Finite element matrix
    u  = data{2};
    Mu = u.dqm(Xmsh);
    
    % Function to be applied
    F  = data{3};
    Fx = F(X);
    if (size(Fx,2) > 1)
        error('mshIntegral.m : unavailable case')
    end
    
    % Integration
    I = Mu' * Wx * Fx;
    
    
%%% FINITE ELEMENT OPERATOR --> \int_{mesh(x)} psi'(x) psi(x) dx 
elseif (length(data) == 3) && isa(data{2},'fem') && isa(data{3},'fem')
    % Mesh integration data
    Xmsh   = data{1};
    [~,Wx] = Xmsh.qud;    
    Wx     = spdiags(Wx,0,length(Wx),length(Wx));
    
    % Finite element matrix
    u  = data{2};
    Mu = u.dqm(Xmsh);
    
    % Finite element matrix
    v  = data{3};
    Mv = v.dqm(Xmsh);
    
    % Integration
    if iscell(Mu) && iscell(Mv)
        I  = sparse(1,1);
        for i = 1:3
            I = I + Mu{i}' * Wx * Mv{i};
        end        
    else
        I = Mu' * Wx * Mv;
    end
    
    
%%% FINITE ELEMENT OPERATOR --> \int_{mesh(x)} psi(x)' f(x) psi(x) dx 
elseif (length(data) == 4) && isa(data{1},'msh') && isa(data{2},'fem')
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx);
    
    % Finite element matrix
    u  = data{2};
    Mu = u.dqm(Xmsh);
    
    % Function to be applied
    F  = data{3};
    Fx = F(X);
    if (size(Fx,2) > 1)
        error('mshIntegral.m : unavailable case')
    else
        Fx = spdiags(Fx,0,Nx,Nx);
    end
    
    % Finite element matrix
    v  = data{4};
    Mv = v.dqm(Xmsh);
    
    % Integration
    if iscell(Mu) && iscell(Mv)
        I  = sparse(1,1);
        for i = 1:3
            I = I + Mu{i}' * Wx * Fx * Mv{i};
        end
    else
        I = Mu' * Wx * Fx * Mv;
    end
    
    
%%% BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(y)} f(x,y) psi(y) dy  
elseif (length(data) == 4) && isnumeric(data{1}) && isa(data{2},'msh')
    % Evaluation points
    X  = data{1};
    Nx = size(X,1);
    
    % Mesh integration data
    Ymsh   = data{2};
    [Y,Wy] = Ymsh.qud;    
    Ny     = size(Y,1);    
    Wy     = spdiags(Wy,0,Ny,Ny); 
    
    % Function to be applied
    F   = data{3};
    Fxy = zeros(Nx,Ny);
    for j = 1:Ny
        Fxy(:,j) = F(X,Y(j,:));
    end
    
    % Finite element matrix
    v  = data{4};
    Mv = v.dqm(Ymsh);
    
    % Integration    
    I = Fxy * Wy * Mv;
    
    
%%% BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x,y) dx  
elseif (length(data) == 4) && isa(data{1},'msh') && isnumeric(data{2})
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;    
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx); 
    
    % Evaluation points
    Y  = data{2};
    Ny = size(Y,1);
    
    % Finite element matrix
    u  = data{3};
    Mu = u.dqm(Xmsh);
   
    % Function to be applied
    F   = data{4};
    Fxy = zeros(Nx,Ny);
    for j = 1:Ny
        Fxy(:,j) = F(X,Y(j,:));
    end
    
    % Integration
    I = Mu' * Wx * Fxy;
    
    
%%% BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy    
elseif (length(data) == 5) && isa(data{1},'msh') && isa(data{2},'msh')
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Nx     = size(X,1);
    Wx     = spdiags(Wx,0,Nx,Nx); 
    
    % Mesh integration data
    Ymsh   = data{2};
    [Y,Wy] = Ymsh.qud;
    Ny     = size(Y,1);
    Wy     = spdiags(Wy,0,Ny,Ny);
    
    % Finite element matrix
    u  = data{3};
    Mu = u.dqm(Xmsh);
    
    % Function applied to quadrature
    F   = data{4};
    Fxy = zeros(Nx,Ny);
    for j = 1:Ny
        Fxy(:,j) = F(X,Y(j,:));
    end
    
    % Finite element matrix
    v  = data{5};
    Mv = v.dqm(Ymsh);
    
    % Integration
    if iscell(Mu) && iscell(Mv)
        I = 0;
        for i = 1:3
            I = I + Mu{i}' * Wx * Fxy * Wy * Mv{i};
        end
    else
        I = Mu' * Wx * Fxy * Wy * Mv;
    end
    
    
%%% H-MATRIX BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(y)} f(x,y) psi(y) dy  
elseif (length(data) == 5) && isnumeric(data{1}) && isa(data{2},'msh')
    % Evaluation points
    X  = data{1};
    Nx = size(X,1);
    
    % Collocation identity matrix
    Mx = speye(Nx,Nx);
    
    % Mesh integration data
    Ymsh   = data{2};
    [Y,Wy] = Ymsh.qud;
    Wy     = spdiags(Wy,0,length(Wy),length(Wy));
    
    % Green kernel
    green = data{3};
    
    % Integrated finite element matrix
    v  = data{4};
    Mv = v.dqm(Ymsh);
    Mv = Wy * Mv;
    
    % Accuracy
    tol = data{5};
   
    % H-Matrix integration
    I = hmx(X,v.dof(Ymsh),Mx,X,green,Y,Mv,tol);


%%% H-MATRIX BOUNDARY ELEMENT INTEGRATION --> \int_{mesh(x)} psi(x)' f(x,y) dx  
elseif (length(data) == 5) && isa(data{1},'msh') && isnumeric(data{2}) 
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Wx     = spdiags(Wx,0,length(Wx),length(Wx));
    
    % Evaluation points
    Y  = data{2};
    Ny = size(Y,1);
    
    % Collocation identity matrix
    My = speye(Ny,Ny);
    
    % Integrated finite element matrix
    u  = data{3};
    Mu = u.dqm(Xmsh);
    Mu = Wx * Mu;
      
    % Green kernel
    green = data{4};
    
    % Accuracy
    tol = data{5};
   
    % H-Matrix integration
    I = hmx(u.dof(Xmsh),Y,Mu,X,green,Y,My,tol);

    
%%% H-MATRIX BOUNDARY ELEMENT OPERATOR  --> \int_{mesh(x)} \int_{mesh(y)} psi(x)' f(x,y) psi(y) dxdy    
elseif (length(data) == 6)
    % Mesh integration data
    Xmsh   = data{1};
    [X,Wx] = Xmsh.qud;
    Wx     = spdiags(Wx,0,length(Wx),length(Wx));
    
    % Mesh integration data
    Ymsh   = data{2};
    [Y,Wy] = Ymsh.qud;
    Wy     = spdiags(Wy,0,length(Wy),length(Wy));
    
    % Finite element matrix
    u  = data{3};
    Mu = u.dqm(Xmsh);
    
    % Green kernel
    green = data{4};
    
    % Finite element
    v  = data{5};
    Mv = v.dqm(Ymsh);
    
    % Integration
    if iscell(Mu) && iscell(Mv)
        for i = 1:3
            Mu{i} = Wx * Mu{i};
            Mv{i} = Wy * Mv{i};
        end
    else
       Mu = Wx * Mu;
       Mv = Wy * Mv;
    end
    
    % Accuracy
    tol = data{6};
   
    % H-Matrix integration
    I = hmx(u.dof(Xmsh),v.dof(Ymsh),Mu,X,green,Y,Mv,tol);
    
    
% UNKNOWN
else
    error('mshIntegral.m : unavailable case')
end

end