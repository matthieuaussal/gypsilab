function Ms = domRegularize(data)
%+========================================================================+
%|                                                                        |
%|              OPENDOM - LIBRARY FOR NUMERICAL INTEGRATION               |
%|           openDom is part of the GYPSILAB toolbox for Matlab           |
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
%|    #    |   FILE       : domRegularize.m                               |
%|    #    |   VERSION    : 0.30                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal & François Alouges            |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 05.09.2017                                    |
%| ( === ) |   SYNOPSIS   : Finite element regularization matrix for      |
%|  `---'  |                singularities with Laplace kernel             |
%+========================================================================+

%%% INPUT ANALYSIS
if (length(data) == 4)
    X     = data{1};
    Ydom  = data{2};
    green = data{3};
    v     = data{4};
else
    Xdom  = data{1};
    Ydom  = data{2};
    u     = data{3};
    green = data{4};
    v     = data{5};
end


%%% INITIALIZATION
% Mesh data from Y
vtx  = Ydom.msh.vtx;
elt  = Ydom.msh.elt;
ctr  = Ydom.msh.ctr;
nrm  = Ydom.msh.nrm;
srf  = Ydom.msh.ndv;
tau  = cell2mat(Ydom.msh.tgt);
nu   = cell2mat(Ydom.msh.nrmEdg);
Nelt = size(elt,1);

% Quadrature data from Y
[Y,Wy,elt2qud] = Ydom.qud;
Yun            = ones(1,size(elt2qud,2));

% Degrees of freedom from Y
[~,elt2dof] = v.dof;
Nbas        = size(elt2dof,2);

% Quadrature data from X
if (length(data) == 5)
    [X,Wx] = Xdom.qud;
end
Nx = size(X,1);

% Rangesearch with max(|edge|)_Y
[Ielt,Relt] = knnsearch(X,ctr,'K',50);                        %%% DEBUG %%%
Mx          = cell(Nelt,1);


%%% RIGHT INTEGRATION WITH REGULARIZATION
for el = 1:Nelt
    % Triangular data for Y
    Sel  = vtx(elt(el,:),:);
    Nel  = nrm(el,:);
    Tel  = reshape(tau(el,:),3,3)';
    NUel = reshape(nu(el,:),3,3)';
    
    % Local size
    edga = Sel(2,:) - Sel(1,:);
    edgb = Sel(3,:) - Sel(1,:);
    edgc = Sel(3,:) - Sel(2,:);
    rMin = max([norm(edga),norm(edgb),norm(edgc)]);           %%% DEBUG %%%
    
    % Quadratures points in interaction
    Iy = elt2qud(el,:);
    Ix = sort(Ielt(el,Relt(el,:)<rMin))';
    
    %     % Graphical representation                          %%% DEBUG %%%
    %     figure(10)
    %     plot(Ydom.msh.sub(el))
    %     hold on
    %     plot3(X(Ix,1),X(Ix,2),X(Ix,3),'*')
    %     axis equal
    %     hold off
    %     pause(1)
    
    % If interactions
    if ~isempty(Ix)
        %%% CORRECTION WITH SEMI-ANALYTIC INTEGRATION
        % Analytical integration
        [Rm1,rRm1,gradRm1] = domSemiAnalyticInt(X(Ix,:),Sel,Nel,Tel,NUel,1e-8);
        %         Rm1(:) = 0; rRm1(:) = 0; gradRm1(:) = 0;    %%% DEBUG %%%
        
        % Vector yg-x
        Xun = ones(length(Ix),1);
        XY1 = Xun * Y(Iy,1)' - X(Ix,1) * Yun;
        XY2 = Xun * Y(Iy,2)' - X(Ix,2) * Yun;
        XY3 = Xun * Y(Iy,3)' - X(Ix,3) * Yun;
        
        % Distance r = |yg-x|
        Rxy             = sqrt(XY1.^2 + XY2.^2 + XY3.^2);
        Rxym1           = 1./Rxy;
        Rxym1(Rxy<1e-6) = 0;
        
        % Int_el(1/|r|) - Sum_g 1/|yg-x|
        Rm1 = Rm1 - Rxym1 * Wy(Iy);
        %         norm(Rm1)                                   %%% DEBUG %%%
        
        % Int_el(r/|r|) - Sum_g (yg-x)/|yg-x|
        rRm1(:,1) = rRm1(:,1) - (XY1 .* Rxym1) * Wy(Iy);
        rRm1(:,2) = rRm1(:,2) - (XY2 .* Rxym1) * Wy(Iy);
        rRm1(:,3) = rRm1(:,3) - (XY3 .* Rxym1) * Wy(Iy);
        %         norm(rRm1)                                  %%% DEBUG %%%
        
        % Int_el(-r/|r|^3) - Sum_g -(yg-x)/|yg-x|^3
        gradRm1(:,1) = gradRm1(:,1) + (XY1 .* Rxym1.^3) * Wy(Iy);
        gradRm1(:,2) = gradRm1(:,2) + (XY2 .* Rxym1.^3) * Wy(Iy);
        gradRm1(:,3) = gradRm1(:,3) + (XY3 .* Rxym1.^3) * Wy(Iy);
        %         norm(gradRm1)                               %%% DEBUG %%%
        
        % Nullify V
        V = [];
        
        %%% FINITE ELEMENT P0
        if strcmp(v.typ,'P0')
            % Correction
            if strcmp(green,'[1/r]') && strcmp(v.opr,'[psi]')
                V = Rm1;
                
            elseif strcmp(green,'[1/r]') && strcmp(v.opr,'n*[psi]')
                V{1} = Rm1 .* Nel(1);
                V{2} = Rm1 .* Nel(2);
                V{3} = Rm1 .* Nel(3);
                
            elseif strcmp(green,'[1/r]') && strcmp(v.opr,'nxgrad[psi]')
                V{1} = zeros(size(Rm1));
                V{2} = zeros(size(Rm1));
                V{3} = zeros(size(Rm1));
                
            elseif strcmp(green,'grady[1/r]') && strcmp(v.opr,'n*[psi]')
                V = gradRm1 * Nel';
                
            else
                error('domRegularize : unavailable case')
            end
            
            
        %%% FINITE ELEMENT P1
        elseif strcmp(v.typ,'P1')
            % For each basis function
            for j = 1:Nbas
                % Next dof
                jp1 = mod(j,3) + 1;
                
                % Height from j
                hj = ((Sel(j,:)-Sel(jp1,:)) * NUel(j,:)');
                
                % Scalar product (x-yk).nuj/hj
                tmp = ( (X(Ix,1)-Sel(jp1,1))*NUel(j,1) + ...
                    (X(Ix,2)-Sel(jp1,2))*NUel(j,2) + ...
                    (X(Ix,3)-Sel(jp1,3))*NUel(j,3) ) ./ hj;
                
                % Correction
                if strcmp(green,'[1/r]') && strcmp(v.opr,'[psi]')
                    V(:,j) = Rm1.*tmp + rRm1*NUel(j,:)'/hj;
                    
                elseif strcmp(green,'[1/r]') && strcmp(v.opr,'n*[psi]')
                    Vx        = Rm1.*tmp + rRm1*NUel(j,:)'/hj;
                    V{1}(:,j) = Vx .* Nel(1);
                    V{2}(:,j) = Vx .* Nel(2);
                    V{3}(:,j) = Vx .* Nel(3);
                    
                elseif strcmp(green,'[1/r]') && strcmp(v.opr,'nxgrad[psi]')
                    NxNUj     = cross(Nel,NUel(j,:));
                    V{1}(:,j) = NxNUj(1)/hj .* Rm1;
                    V{2}(:,j) = NxNUj(2)/hj .* Rm1;
                    V{3}(:,j) = NxNUj(3)/hj .* Rm1;
                    
                elseif strcmp(green,'grady[1/r]') && strcmp(v.opr,'n*[psi]')
                    V(:,j) = tmp .* (gradRm1 * Nel');
                    
                else
                    error('domRegularize : unavailable case')
                end
            end
            
            
        %%% FINITE ELEMENT RWG
        elseif strcmp(v.typ,'RWG')
            % For each basis function
            for j = 1:Nbas
                % Flux inside element
                jp1  = mod(j,3) + 1;
                jp2  = mod(jp1,3) + 1;
                sgnS = (2*(elt(el,jp1) < elt(el,jp2))-1)./(2*srf(el));
                
                % Vertex to quadrature point
                XmS = X(Ix,:) - Xun * Sel(j,:);
                
                % Correction
                if strcmp(green,'[1/r]') && strcmp(v.opr,'[psi]')
                    V{1}(:,j) = sgnS .* (XmS(:,1).*Rm1 + rRm1(:,1));
                    V{2}(:,j) = sgnS .* (XmS(:,2).*Rm1 + rRm1(:,2));
                    V{3}(:,j) = sgnS .* (XmS(:,3).*Rm1 + rRm1(:,3));
                    
                elseif strcmp(green,'[1/r]') && strcmp(v.opr,'div[psi]')
                    V(:,j) = 2*sgnS .* Rm1;
                    
                elseif strcmp(green,'grady[1/r]') && strcmp(v.opr,'[psi]')
                    gradRm1xXmS = cross(gradRm1,XmS);
                    V{1}(:,j)   = sgnS .* gradRm1xXmS(:,1);
                    V{2}(:,j)   = sgnS .* gradRm1xXmS(:,2);
                    V{3}(:,j)   = sgnS .* gradRm1xXmS(:,3);
                    
                else
                    error('domRegularize : unavailable case')
                end
                
            end
            
        else
            error('domRegularize : unavailable case')
        end
        
        % Matrix-Vector product
        I = Ix * ones(1,Nbas);
        J = ones(length(Ix),1) * elt2dof(el,:);
        if iscell(V)
            Mx{el} = [I(:) J(:) V{1}(:) V{2}(:) V{3}(:)];
        else
            Mx{el} = [I(:) J(:) V(:)];
        end
    end
end


%%% LEFT INTEGRATION AND RIGHT CORRECTIONS
% Right unknowns
[Xunk,Punk] = v.unk;

% Data preparation
Ndof = size(v.dof,1);
Mx   = cell2mat(Mx);
if (length(data) == 4)
    Mu = speye(Nx,Nx);
    Mw = speye(Nx,Nx);
    Ms = sparse(Nx,Ndof);
else
    Mu = u.dqm(Xdom);
    Mw = spdiags(Wx,0,length(Wx),length(Wx));
    Ms = sparse(size(Xunk,1),Ndof);
end

% Integration
if iscell(Mu)
    for i = 1:length(Mu)
        Ms = Ms + Mu{i}' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,2+i),Nx,Ndof);
    end
else
    Ms = Mu' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,3),Nx,Ndof);
end

% Right reduction
Ms = Ms * Punk;
end
