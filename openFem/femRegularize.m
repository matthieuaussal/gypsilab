function Ms = femRegularize(varargin)
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
%|              François Alouges - CMAP, Ecole Polytechnique              |
%| Creation   : 14.03.17                                                  |
%| Last modif : 21.06.17                                                  |
%| Synopsis   : Finite element regularization matrix for laplace kernel   |
%+========================================================================+

%%% INPUT ANALYSIS
if (nargin == 4)
    X     = varargin{1};
    Ymsh  = varargin{2};
    green = varargin{3};
    v     = varargin{4};
else
    Xmsh  = varargin{1};
    Ymsh  = varargin{2};
    u     = varargin{3};
    green = varargin{4};
    v     = varargin{5};
end


%%% INITIALIZATION
% Mesh data from Y
vtx  = Ymsh.vtx;
elt  = Ymsh.elt;
ctr  = Ymsh.ctr;
nrm  = Ymsh.nrm;
tau  = cell2mat(Ymsh.tgt);
nu   = cell2mat(Ymsh.nrmEdg);
Nelt = size(elt,1);

% Quadrature data from Y
[Y,Wy,elt2qud] = Ymsh.qud;
Yun            = ones(1,size(elt2qud,2));

% Degrees of freedom from Y
[~,elt2dof] = v.dof(Ymsh);
Nbas        = size(elt2dof,2);

% Quadrature data from X
if (nargin == 5)
    [X,Wx] = Xmsh.qud;
end
Nx = size(X,1);

% Rangesearch with max(|edge|)_Y
[Ielt,Relt] = knnsearch(X,ctr,'K',50);
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
    rMin = max([norm(edga),norm(edgb),norm(edgc)]);
    
    % Quadratures points in interaction
    Iy = elt2qud(el,:);
    Ix = sort(Ielt(el,Relt(el,:)<rMin))';

%     % Graphical representation
%     figure(10)
%     plot(Ymsh.sub(el))
%     hold on
%     plot3(X(Ix,1),X(Ix,2),X(Ix,3),'*')
%     axis equal
%     hold off
%     pause(1)

    % If interactions
    if ~isempty(Ix)
        %%% CORRECTION WITH SEMI-ANALYTIC INTEGRATION        
        % Analytical integration
        [Rm1,rRm1,gradRm1] = mshSemiAnalyticInt(X(Ix,:),Sel,Nel,Tel,NUel,1e-8);
%         Rm1(:) = 0; rRm1(:) = 0; gradRm1(:) = 0;
        
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
%         norm(Rm1)
        
        % Int_el(r/|r|) - Sum_g (yg-x)/|yg-x|
        rRm1(:,1) = rRm1(:,1) - (XY1 .* Rxym1) * Wy(Iy);
        rRm1(:,2) = rRm1(:,2) - (XY2 .* Rxym1) * Wy(Iy);
        rRm1(:,3) = rRm1(:,3) - (XY3 .* Rxym1) * Wy(Iy);
%         norm(rRm1)
        
        % Int_el(-r/|r|^3) - Sum_g -(yg-x)/|yg-x|^3
        gradRm1(:,1) = gradRm1(:,1) + (XY1 .* Rxym1.^3) * Wy(Iy);
        gradRm1(:,2) = gradRm1(:,2) + (XY2 .* Rxym1.^3) * Wy(Iy);
        gradRm1(:,3) = gradRm1(:,3) + (XY3 .* Rxym1.^3) * Wy(Iy);
%         norm(gradRm1)
        
       
        %%% FINITE ELEMENT P0
        if strcmp(v.typ,'P0')
            % Quadrature
            if strcmp(green,'[1/r]')
                Vx = Rm1;
                
            elseif strcmp(green(1:end-1),'gradx[1/r]')
                d  = str2double(green(end));
                Vx = - gradRm1(:,d);
                
            elseif strcmp(green(1:end-1),'grady[1/r]')
                d  = str2double(green(end));
                Vx = gradRm1(:,d);
                
            else
                error('femRegularize : unavailable case')
            end
            
            % Finite element
            if strcmp(v.opr,'[psi]')
                V = Vx;
            
            elseif strcmp(v.opr,'n*[psi]')
                V{1} = Vx .* Nel(1);
                V{2} = Vx .* Nel(2);
                V{3} = Vx .* Nel(3);
                
            elseif strcmp(v.opr,'nxgrad[psi]')
                V{1} = zeros(size(Vx));
                V{2} = zeros(size(Vx));
                V{3} = zeros(size(Vx));
                
            elseif strcmp(v.opr(1:end-1),'n*[psi]')
                d = str2double(v.opr(end));
                V = Vx .* Nel(d);
                
            elseif strcmp(v.opr(1:end-1),'nxgrad[psi]')
                V = zeros(size(Vx));
                
            else
                error('femRegularize : unavailable case')
            end   
            
            
        %%% FINITE ELEMENT P1
        elseif strcmp(v.typ,'P1')
            % Initialization
            Vx = zeros(length(Ix),Nbas);
            if strcmp(v.opr,'n*[psi]') || strcmp(v.opr,'nxgrad[psi]')
                V = cell(1,3);
            else
                V = zeros(length(Ix),Nbas);
            end
            
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
                
                % Quadrature
                if strcmp(green,'[1/r]')
                    Vx(:,j) = Rm1.*tmp + rRm1*NUel(j,:)'/hj;
                    
                elseif strcmp(green(1:end-1),'gradx[1/r]')
                    d       = str2double(green(end));
                    Vx(:,j) = - tmp .* gradRm1(:,d);
                    
                elseif strcmp(green(1:end-1),'grady[1/r]')
                    d       = str2double(green(end));
                    Vx(:,j) = tmp .* gradRm1(:,d);
                    
                else
                    error('femRegularize : unavailable case')
                end
                
                % Finite element
                if strcmp(v.opr,'[psi]')
                    V = Vx;
                    
                elseif strcmp(v.opr,'n*[psi]')
                    V{1}(:,j) = Vx(:,j) .* Nel(1);
                    V{2}(:,j) = Vx(:,j) .* Nel(2);
                    V{3}(:,j) = Vx(:,j) .* Nel(3);
                    
                elseif strcmp(v.opr,'nxgrad[psi]')
                    V{1}(:,j) = 1/hj*(Nel(2)*NUel(j,3)-Nel(3)*NUel(j,2)) .* Rm1;
                    V{2}(:,j) = 1/hj*(Nel(3)*NUel(j,1)-Nel(1)*NUel(j,3)) .* Rm1;
                    V{3}(:,j) = 1/hj*(Nel(1)*NUel(j,2)-Nel(2)*NUel(j,1)) .* Rm1;
                    
                elseif strcmp(v.opr(1:end-1),'n*[psi]')
                    d      = str2double(v.opr(end));
                    V(:,j) = Vx(:,j) .* Nel(d);
                    
                elseif strcmp(v.opr(1:end-1),'nxgrad[psi]')
                    d      = str2double(v.opr(end));
                    dp1    = mod(d,3) + 1;
                    dp2    = mod(dp1,3) + 1;
                    V(:,j) = 1/hj*(Nel(dp1)*NUel(j,dp2)-Nel(dp2)*NUel(j,dp1)) .* Rm1;
                    
                else
                    error('femRegularize : unavailable case')
                end
                
            end
        else
            error('femRegularize : unavailable case')
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


%%% LEFT INTEGRATION
% Data preparation
Ndof = size(v.dof(Ymsh),1);
Mx   = cell2mat(Mx);
if (nargin == 4)
    Mu = speye(Nx,Nx);
    Mw = speye(Nx,Nx);
    Ms = sparse(Nx,Ndof);
else
    Mu = u.dqm(Xmsh);
    Mw = spdiags(Wx,0,length(Wx),length(Wx));
    Ms = sparse(size(u.dof(Xmsh),1),Ndof);
end

% Integration
if iscell(Mu)
    for i = 1:length(Mu)
        Ms = Ms + Mu{i}' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,2+i),Nx,Ndof); 
    end
else
    Ms = Mu' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,3),Nx,Ndof);
end

end
