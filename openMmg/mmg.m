classdef mmg < handle
% Copyright (c) 2018-2019
% Matthieu Aussal, CMAP, Ecole Polytechnique  
% Algiane Froehly, CARDAMOME, INRIA-SOFT 
% LGPL Lesser General Public License v3.0. 
% Remeshing using Mmg tools : https://www.mmgtools.org             

properties 
    % EDGE 
    Aniso    = [];   % AUTOMATIC ANISOTROPY
    Optim    = [];   % AUTOMATIC MESH IMPROVEMENT
    Hgrad    = [];   % MAXIMAL RATIO BETEEN 2 ADJACENT EDGES
    Hgradreq = [];   % RATIO BETWEEN REQUIRED ENTITIES AND NEIGHBOURG 
    Hmax     = [];   % MAXIMAL EDGE LENGTH
    Hmin     = [];   % MINIMAL EDGE LENGTH
    Hsiz     = [];   % CONSTANT EDGE LENGTH
    Hausd    = [];   % HAUSSDORF DISTANCE
    
    % Adaptation
    Nofem    = [];   % REMESHING WITHOUT FINITE ELEMENT ENFORCEMENT 
    Noinsert = [];   % CONSTANT NODES NUMBER
    Nomove   = [];   % NO POINT RELOCATION 
    Nosurf   = [];   % NO SURFACE MODIFICATION  
    Noswap   = [];   % NO EDGE SWAPPING
    Angle    = [];   % SHARP ANGLE DETECTION (°) 
    
    % General
    Memory  = [];    % MAXIMUM MEMORY USED (M0)
    Verbose = [];    % VERBOSITY
    Options = [];    % PRINT OPTIONS DEFINITION
    
    % Mesh and size map
    Mesh    = [];    % MESH OBJECT FROM GYPSILAB (vtx, elt, col)
    Map     = [];    % SIZE MAP AT VERTICES
end

methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function obj = mmg(varargin)
        % Initialize
        if (nargin==0)
            
        % Input is mesh    
        elseif (nargin==1)
            obj.Mesh = varargin{1};
            
        % Input is mesh and Haussdorf distance
        elseif (nargin==2)
            obj.Mesh  = varargin{1};
            obj.Hausd = ['-hausd ',num2str(varargin{2})];
            
        else
            error('mmg.m : unavailable case')
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACCESSORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function aniso(mmg)
        mmg.Aniso = '-A';
    end
    
    function optim(mmg)
        mmg.Optim = '-optim';
    end
    
    function hgrad(mmg,val)
        mmg.Hgrad = ['-hgrad ',num2str(val)];
    end
    
    function hgradreq(mmg,val)
        mmg.Hgradreq = ['-hgradreq ',num2str(val)];
    end
    
    function hmax(mmg,val)
        mmg.Hmax = ['-hmax ',num2str(val)];
    end
    
    function hmin(mmg,val)
        mmg.Hmin = ['-hmin ',num2str(val)];
    end
    
    function hsiz(mmg,val)
        mmg.Hsiz  = ['-hsiz ',num2str(val)];
    end
    
    function hausd(mmg,val)
        mmg.Hausd  = ['-hausd ',num2str(val)];
    end
    
    function nofem(mmg)
        mmg.Nofem = '-nofem';
    end
    
    function noinsert(mmg)
        mmg.Noinsert = '-noinsert';
    end
    
    function nomove(mmg)
        mmg.Nomove = '-nomove';
    end
    
    function nosurf(mmg)
        mmg.Nosurf = '-nosurf';
    end
    
    function noswap(mmg)
        mmg.Noswap = '-noswap';
    end
    
    function angle(mmg,val)
        mmg.Angle = ['-ar ',num2str(val)];
    end
    
    function memory(mmg,val)
        mmg.Memory = ['-m ',num2str(val)];
    end
    
    function verbose(mmg,val)
        mmg.Verbose = ['-v ',num2str(val)];
    end
    
    function options(mmg)
        mmg.Options = '-h ';
    end
        
    function mesh(mmg,msh)
        mmg.Mesh = msh;
    end
    
    function map(mmg,map)
        mmg.Map = map;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [mesh,val] = run(mmg)
        % Current and mmg directory
        here  = pwd;
        there = which('mmg.m');
        there = there(1:end-6);
        
        % Move to mmg directory
        cd(there)

        % Check mesh (3D nodes, triangle or tetra, colours are int)
        mesh = mmg.Mesh;
        if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) < 3) || (size(mesh.elt,2) > 4) 
            error('mmg.m : unavailable case, please use 3D vertices.')
        end
        
        % Check colours are int
        if (norm(mesh.col-floor(mesh.col),'inf') > 1e-12) || (min(mesh.col)<0)
            error('mmg.m : only integer values for colours.')
        end
        
        % Write original mesh in .msh format
        if isempty(mmg.Map)
            mshWriteMsh('original.msh',mesh);
        else
            mshWriteMsh('original.msh',mesh,mmg.Map);
        end
            
        % Operating system
        if ismac
            os = './mac_';
        elseif ispc
            os = 'win_';
        elseif isunix
            os = './uni_';
        else
            disp('mmg.m : unavailable case')
        end
        
        % Element form
        if (size(mesh.elt,2) == 3)
            bin = 'mmgs_O3';
        elseif (size(mesh.elt,2) == 4)
            bin = 'mmg3d_O3';
        end
        
        % Add extension (for windows only)
        if ispc
            bin = [bin,'.exe '];
        else
            bin = [bin,' '];
        end
        
        % Define temporary mesh files
        file = '-in original.msh -out refined.msh ';
        
        % Convert input to command
        field = fieldnames(mmg);
        opt   = '';
        for i = 1:length(field)
            cmd = getfield(mmg,field{i});
            if ~isempty(cmd) && ~strcmp(field{i},'Mesh') && ~strcmp(field{i},'Map')
                opt = [opt,cmd,' '];
            end                
        end
        
        % Execute mmg binaries
        command = [os,bin,file,opt];
        system(command);
        
        % Read refined mesh
        [vtx,elt,col,val] = mshReadMsh('refined.msh');
        mesh = msh(vtx,elt,col);
        
        % Clean meshes
        delete('original.msh')
        delete('refined.msh')
        
        % Back to current directory
        cd(here)
    end
end
end
