classdef msh
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
%| Synopsis   : Mesh class definition                                     |
%+========================================================================+

properties
    vtx = [];      % VERTEX COORDINATES (3 dimensions)
    elt = [];      % ELEMENTS LIST (particles, edges, triangles or tetrahedron)
    gss = 1;       % NUMBER OF INTEGRATION POINTS
end

methods
    % CONSTRUCTOR
    function mesh = msh(varargin)
        % Read file
        if (length(varargin) == 1)
            mesh = msh;
            ext  = varargin{1}(end-2:end);
            if strcmp(ext,'ply')
                [mesh.vtx,mesh.elt] = mshReadPly(varargin{1});
            elseif strcmp(ext,'vtk')
                [mesh.vtx,mesh.elt] = mshReadVtk(varargin{1});
            else
                error('msh.m : unavailable case')
            end
        
        % Use input data    
        elseif (length(varargin) == 2)
            mesh.vtx = varargin{1};
            mesh.elt = varargin{2};
        end
        
        % Clean mesh
        mesh = mshClean(mesh);
    end
    
    % FREE BOUNDARY
    function bnd2vtx = bnd(mesh)
        bnd2vtx = mshBoundary(mesh);
    end
    
    % CENTER
    function X = ctr(mesh)
        X = mshCenter(mesh);
    end
    
    % EDGES
    function [edg2vtx,elt2edg] = edg(mesh)
        [edg2vtx,elt2edg] = mshEdge(mesh);
    end
    
    % EDGES LENGTH
    function lgt = edgLgt(mesh)
        edg2vtx = edg(mesh);
        lgt = mesh.vtx(edg2vtx(:,2),:) - mesh.vtx(edg2vtx(:,1),:);
        lgt = sqrt(sum(lgt.^2,2));
        lgt = [min(lgt) max(lgt)];
    end
    
    % FACES
    function [fce2vtx,elt2fce] = fce(mesh)
        [fce2vtx,elt2fce] = mshFace(mesh);
    end
    
    % NUMERICAL INTEGRATION
    function I = integral(varargin)
        I = mshIntegral(varargin);
    end  
    
    % MESH REPRESENTATION
    function plot(varargin)
        mshPlot(varargin);
    end
    
    % GAUSSIAN QUADRATURE
    function [X,W,elt2qud] = qud(mesh)
        [X,W,elt2qud] = mshQuadrature(mesh);
    end
        
    % SUBMESHING
    function mesh = sub(mesh,Ielt)
        mesh.elt = mesh.elt(Ielt,:);
        mesh     = mshClean(mesh);
    end
        
    % SURFACES
    function S = ndv(mesh)
        S = mshNdvolume(mesh);
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%% TRIANGULAR MESH ONLY %%%%%%%%%%%%%%%%%%%%%%%%%    
    % TANGENTS
    function T = tgt(mesh)
        if (size(mesh.elt,2) == 3)
            T = cell(1,3);
            for i = 1:3
                ip1  = mod(i,3)+1;
                ip2  = mod(ip1,3)+1;
                A    = mesh.vtx(mesh.elt(:,ip1),:);
                B    = mesh.vtx(mesh.elt(:,ip2),:);
                T{i} = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
            end
        else
            error('msh.m : unavailable case')
        end
    end

    % NORMALS
    function N = nrm(mesh)
        T = mesh.tgt;
        N = cross(T{1},T{2});
        N = N ./ (sqrt(sum(N.^2,2)) * [1 1 1]);
    end
    
    % EDGES NORMALS
    function Nu = nrmEdg(mesh)
        Nu = cell(1,3);
        for i = 1:3
            Nu{i} = cross(mesh.tgt{i},mesh.nrm);
        end
    end
    
    % GAUSSIAN NORMALES
    function N = nrmQud(mesh)
        x = mshReference(mesh);
        N = zeros(size(x,1)*size(mesh.elt,1),size(mesh.vtx,2));
        for j = 1:size(x,1)
            idx      = (j:size(x,1):size(x,1) * size(mesh.elt,1))';
            N(idx,:) = mesh.nrm;
        end
    end

end
end
