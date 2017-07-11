function M = femLagrangePn(fe,mesh)
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
%| Synopsis   : Finite element matrices                                   |
%+========================================================================+

%%% FINITE ELEMENT MATRIX AND GRADIENT
if strcmp(fe.opr,'[psi]') || strcmp(fe.opr,'grad[psi]')
    % Gaussian quadrature
    x                = mshReference(mesh);
    [Xqud,~,elt2qud] = mesh.qud;
    
    % Degrees of freedom
    [Xdof,elt2dof] = fe.dof(mesh);
    
    % Dimensions
    Nelt = size(mesh.elt,1);
    Nqud = size(Xqud,1);
    Ngss = size(elt2qud,2);
    Ndof = size(Xdof,1);
    Nbas = size(elt2dof,2);
    
    % Basis fct constant per elements
    if strcmp(fe.typ,'P0')
        % Edge mesh
        if (size(mesh.elt,2) == 2)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);

        % Triangular mesh
        elseif (size(mesh.elt,2) == 3)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);
            dyF = zeros(1,Ngss);
            
        % Tetrahedron mesh
        elseif (size(mesh.elt,2) == 4)
            F   = ones(1,Ngss);
            dxF = zeros(1,Ngss);
            dyF = zeros(1,Ngss);
            dzF = zeros(1,Ngss);
        
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    % Basis fct piecewise linear per element
    elseif strcmp(fe.typ,'P1')
        % Initialization
        F   = zeros(size(mesh.elt,2),Ngss);
        dxF = zeros(size(mesh.elt,2),Ngss);
        dyF = zeros(size(mesh.elt,2),Ngss);
        dzF = zeros(size(mesh.elt,2),Ngss) ;       
        
        % Edge mesh
        if (size(mesh.elt,2) == 2)
            F(1,:) = 1 - x(:,1);
            F(2,:) = x(:,1);
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;

        % Triangular mesh
        elseif (size(mesh.elt,2) == 3)
            F(1,:) = 1 - x(:,1) - x(:,2);
            F(2,:) = x(:,1);
            F(3,:) = x(:,2);            
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;
            
            dyF(1,:) = - 1;
            dyF(3,:) = 1;
            
        % Tetrahedron mesh
        elseif (size(mesh.elt,2) == 4)
            F(1,:) = 1 - x(:,1) - x(:,2) - x(:,3);
            F(2,:) = x(:,1);
            F(3,:) = x(:,2);
            F(4,:) = x(:,3);
            
            dxF(1,:) = - 1;
            dxF(2,:) = 1;
                        
            dyF(1,:) = - 1;
            dyF(3,:) = 1;
            
            dzF(1,:) = - 1;
            dzF(4,:) = 1;
            
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    % Lagrange order 2, piecewise quadratic by element     
    elseif strcmp(fe.typ,'P2')
        % Edge mesh
        if (size(mesh.elt,2) == 2)
            % Initialization
            F   = zeros(3,Ngss);
            dxF = zeros(3,Ngss);
            dyF = zeros(3,Ngss);
            dzF = zeros(3,Ngss);
            X = x(:,1);
            
            F(1,:) = (1 - X).*(1 - 2*X);
            F(2,:) = X.*(2*X - 1);
            F(3,:) = X.*(1 - X);
            
            dxF(1,:) = 4*X - 3;
            dxF(2,:) = 4*X - 1;
            dxF(3,:) = 2*X - 1;
            
        % Triangular mesh
        elseif (size(mesh.elt,2) == 3)
            % Initialization
            F   = zeros(6,Ngss);
            dxF = zeros(6,Ngss);
            dyF = zeros(6,Ngss);
            dzF = zeros(6,Ngss);
            X = x(:,1);
            Y = x(:,2);
            
            F(1,:) = (1 - X - Y).*(1 - 2*X - 2*Y);
            F(2,:) = X.*(2*X - 1);
            F(3,:) = Y.*(2*Y - 1);
            F(4,:) = 4*X.*(1 - X - Y);
            F(5,:) = 4*X.*Y;
            F(6,:) = 4*Y.*(1 - X - Y);
            
            dxF(1,:) = -3 + 4*(X + Y);
            dxF(2,:) = 4*X - 1;
            dxF(4,:) = 4*(1 - 2*X - Y);
            dxF(5,:) = 4*Y;
            dxF(6,:) = -4*Y;
            
            dyF(1,:) = -3 + 4*(X + Y);
            dyF(3,:) = 4*Y - 1;
            dyF(4,:) = -4*X;
            dyF(5,:) = 4*X;
            dyF(6,:) = 4*(1 - X - 2*Y);
            
        % Tetrahedron mesh
        elseif (size(mesh.elt,2) == 4)
            % Initialization
            F   = zeros(10,Ngss);
            dxF = zeros(10,Ngss);
            dyF = zeros(10,Ngss);
            dzF = zeros(10,Ngss);
            X = x(:,1);
            Y = x(:,2);
            Z = x(:,3);
            
            F(1,:)  = (1 - X - Y - Z).*(1 - 2*X - 2*Y - 2*Z);
            F(2,:)  = X.*(2*X - 1);
            F(3,:)  = Y.*(2*Y - 1);
            F(4,:)  = Z.*(2*Z - 1);
            F(5,:)  = 4*X.*(1 - X - Y - Z);
            F(6,:)  = 4*X.*Y;
            F(7,:)  = 4*Y.*(1 - X - Y - Z);
            F(8,:)  = 4*Z.*(1 - X - Y - Z);
            F(9,:)  = 4*X.*Z;
            F(10,:) = 4*Y.*Z;
             
            dxF(1,:) = -3 + 4*(X + Y + Z);
            dxF(2,:) = 4*X - 1;
            dxF(5,:) = 4*(1 - 2*X - Y - Z);
            dxF(6,:) = 4*Y;
            dxF(7,:) = -4*Y;
            dxF(8,:) = -4*Z;
            dxF(9,:) = 4*Z;
            
            dzF(1,:)  = -3 + 4*(X + Y + Z);
            dzF(4,:)  = 4*Z - 1;
            dzF(5,:)  = -4*X;
            dzF(7,:)  = -4*Y;
            dzF(8,:)  = 4*(1 - X - Y - 2*Z);
            dzF(9,:)  = 4*X;
            dzF(10,:) = 4*Y;
            
            dyF(1,:)  = -3 + 4*(X + Y + Z);
            dyF(3,:)  = 4*Y - 1;
            dyF(5,:)  = -4*X;
            dyF(6,:)  = 4*X;
            dyF(7,:)  = 4*(1 - X - 2*Y - Z);
            dyF(8,:)  = -4*Z;
            dyF(10,:) = 4*Z;

            dzF(1,:)  = -3 + 4*(X + Y + Z);
            dzF(4,:)  = 4*Z - 1;
            dzF(5,:)  = -4*X;
            dzF(7,:)  = -4*Y;
            dzF(8,:)  = 4*(1 - X - Y - 2*Z);
            dzF(9,:)  = 4*X;
            dzF(10,:) = 4*Y;
            
        else
            error('femLagrangePn.m : unavailable case')
        end
        
    else
        error('femLagrangePn.m: unavailable case');
    end
    
    % Numbering dof to quadrature
    idx = zeros(Nbas*Nqud,1);
    jdx = zeros(Nbas*Nqud,1);
    for i = 1:Nbas
        for j = 1:Ngss
            ind = j+(i-1)*Ngss : Nbas*Ngss : Nbas*Nqud;
            idx(ind) = elt2qud(:,j);
            jdx(ind) = elt2dof(:,i);
        end
    end
    
    % Finite element matrix
    if strcmp(fe.opr,'[psi]')
        bas = reshape(F',[1,Ngss*Nbas]);
        bas = kron(ones(1,Nelt),bas);
        M   = sparse(idx,jdx,bas,Nqud,Ndof);
        
    % Gradient of Finite element matrix
    elseif strcmp(fe.opr,'grad[psi]')
        % Edge mesh
        if (size(mesh.elt,2) == 2)
            notYet
            
        % Triangular mesh
        elseif (size(mesh.elt,2) == 3)
            % Vector basis
            E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
            
            % Gramm matrix coefficients
            a = E1(:,1).^2 + E1(:,2).^2 + E1(:,3).^2;
            b = E1(:,1).*E2(:,1) + E1(:,2).*E2(:,2) + E1(:,3).*E2(:,3);
            c = b;
            d = E2(:,1).^2 + E2(:,2).^2 + E2(:,3).^2;
            
            % Determinant
            detG = a.*d - b.*c;
            
            % Inverse by co-factor
            Dx1 = d  ./ detG;
            Dx2 = -b ./ detG;
            Dy1 = -c ./ detG;
            Dy2 = a  ./ detG;
            
            % Gradient projection to integration points
            dbas = cell(1,3);
            for n = 1:3
                dbas{n} = zeros(Nbas,Nqud);
                for d = 1:Nbas
                    DCVx         = Dx1.*E1(:,n) + Dx2.*E2(:,n);
                    DCVy         = Dy1.*E1(:,n) + Dy2.*E2(:,n);
                    dbas{n}(d,:) = kron(DCVx',dxF(d,:)) + kron(DCVy',dyF(d,:));
                end
            end
            
        % Tetrahedron mesh
        elseif (size(mesh.elt,2) == 4)            
            % Vector basis
            E1 = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            E2 = mesh.vtx(mesh.elt(:,3),:) - mesh.vtx(mesh.elt(:,1),:);
            E3 = mesh.vtx(mesh.elt(:,4),:) - mesh.vtx(mesh.elt(:,1),:);
            
            % Gramm matrix
            a = E1(:,1).^2 + E1(:,2).^2 + E1(:,3).^2;
            b = E1(:,1).*E2(:,1) + E1(:,2).*E2(:,2) + E1(:,3).*E2(:,3);
            c = E1(:,1).*E3(:,1) + E1(:,2).*E3(:,2) + E1(:,3).*E3(:,3);
            d = b;
            e = E2(:,1).^2 + E2(:,2).^2 + E2(:,3).^2;
            f = E2(:,1).*E3(:,1) + E2(:,2).*E3(:,2) + E2(:,3).*E3(:,3);
            g = c;
            h = f;
            i = E3(:,1).^2 + E3(:,2).^2 + E3(:,3).^2;
            
            % Determinant (Sarrus rules)
            detG = a.*e.*i + b.*f.*g + c.*d.*h - c.*e.*g - f.*h.*a - i.*b.*d;
            
            % Inverse by co-factor
            Dx1 = (e.*i - f.*h) ./ detG;
            Dx2 = (c.*h - b.*i) ./ detG;
            Dx3 = (b.*f - c.*e) ./ detG;
            Dy1 = (f.*g - d.*i) ./ detG;
            Dy2 = (a.*i - c.*g) ./ detG;
            Dy3 = (c.*d - a.*f) ./ detG;
            Dz1 = (d.*h - e.*g) ./ detG;
            Dz2 = (b.*g - a.*h) ./ detG;
            Dz3 = (a.*e - b.*d) ./ detG;
            
            % Gradient projection to integration points
            dbas = cell(1,3);
            for n = 1:3
                dbas{n} = zeros(Nbas,Nqud);
                for d = 1:Nbas
                    DCVx         = Dx1.*E1(:,n) + Dx2.*E2(:,n) + Dx3.*E3(:,n);
                    DCVy         = Dy1.*E1(:,n) + Dy2.*E2(:,n) + Dy3.*E3(:,n);
                    DCVz         = Dz1.*E1(:,n) + Dz2.*E2(:,n) + Dz3.*E3(:,n);
                    dbas{n}(d,:) = kron(DCVx',dxF(d,:)) ...
                        + kron(DCVy',dyF(d,:)) ...
                        + kron(DCVz',dzF(d,:));
                end
            end
        end
        
        % Dof to quadrature matrix
        M = cell(1,3);
        for n = 1:3
            % Basis functions values
            val = zeros(Nbas*Nqud,1);
            for i = 1:Nbas
                for j = 1:Ngss
                    indl      = j+(i-1)*Ngss : Nbas*Ngss : Nbas*Nqud;
                    indr      = j:Ngss:Nqud;
                    val(indl) = dbas{n}(i,indr);
                end
            end
            
            % Sparse format
            M{n} = sparse(idx,jdx,val,Nqud,Ndof);
        end
    end

    
%%% NORMALS * DQM
elseif strcmp(fe.opr,'n*[psi]')
    % Finite element
    fe.opr = '[psi]';
    dqm    = fe.dqm(mesh);
    
    % Normals
    N = mesh.nrmQud;
    
    % Dot product
    m    = size(N,1);
    M{1} = spdiags(N(:,1),0,m,m) * dqm;
    M{2} = spdiags(N(:,2),0,m,m) * dqm;
    M{3} = spdiags(N(:,3),0,m,m) * dqm;

    
%%%% NORMAL x DQM
elseif strcmp(fe.opr,'nxgrad[psi]')
    % Finite elements
    fe.opr = 'grad[psi]';
    dqm    = fe.dqm(mesh);
    
    % Normals
    nrm  = mesh.nrmQud;
    m    = size(nrm,1);
    N{1} = spdiags(nrm(:,1),0,m,m);
    N{2} = spdiags(nrm(:,2),0,m,m);
    N{3} = spdiags(nrm(:,3),0,m,m);
    
    % Cross product
    M = cell(1,3);
    for i = 1:3
        ip1  = mod(i,3) + 1;
        ip2  = mod(ip1,3) + 1;
        M{i} = N{ip1} * dqm{ip2} - N{ip2} * dqm{ip1};
    end
    
    
%%% GRADIENT (j)
elseif strcmp(fe.opr(1:end-1),'grad[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite elements
    fe.opr = 'grad[psi]';
    dqm    = fe.dqm(mesh);
    
    % Gradient (j)
    M  = dqm{j};

        
%%% NORMALS * DQM (j)
elseif strcmp(fe.opr(1:end-1),'n*[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite element
    fe.opr = '[psi]';
    dqm    = fe.dqm(mesh);
    
    % Normal (component j)
    N = mesh.nrmQud;
    
    % Dot product (j)
    m = size(N,1);
    M = spdiags(N(:,j),0,m,m) * dqm;
    
    
%%%% NORMAL x DQM (j)
elseif strcmp(fe.opr(1:end-1),'nxgrad[psi]')
    % Component
    j = str2double(fe.opr(end));
    
    % Finite element
    fe.opr = 'grad[psi]';
    dqm    = fe.dqm(mesh);
    
    % Normals
    N = mesh.nrmQud;
    m = size(N,1);
    
    % Cross component 
    jp1 = mod(j,3) + 1;
    jp2 = mod(jp1,3) + 1;

    % Cross product
    M = spdiags(N(:,jp1),0,m,m) * dqm{jp2} - ...
        spdiags(N(:,jp2),0,m,m) * dqm{jp1};
    
else
    error('femLagrangePn.m : unavailable case')
end
end
