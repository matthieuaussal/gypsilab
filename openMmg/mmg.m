function [mesh,val] = mmg(varargin)
% Copyright (c) 2018-2019, Matthieu Aussal, Ecole Polytechnique       
% LGPL Lesser General Public License v3.0. 
% Remeshing using level-set from Mmg tools : https://www.mmgtools.org           

% Current and mmg directory
here  = pwd;
there = which('mmg.m');
there = there(1:end-6);

% Move to mmg directory
cd(there)

% Check mesh (3D nodes, triangle or tetra) 
mesh = varargin{1};
if (size(mesh.vtx,2) ~= 3) || (size(mesh.elt,2) < 3) || (size(mesh.elt,2) > 4)
    error('mmg.m : unavailable case')
end

% Centering mesh to origin
ctr      = mean(mesh.vtx,1);
mesh.vtx = mesh.vtx - ones(size(mesh.vtx,1),1)*ctr;

% Normalize mesh (unitary for each direction)
nrm      = mean( sqrt(sum(mesh.vtx.^2,2)) );
mesh.vtx = (1/nrm) .* mesh.vtx;

% Configuration values
if (nargin == 1)   
    stp  = mesh.stp;
    hsiz = stp(3);
elseif (nargin == 2)
    hnod = (1/nrm) .* varargin{2};
elseif (nargin == 3)
    hmin = (1/nrm) .* varargin{2};
    hmax = (1/nrm) .* varargin{3};  
else
    disp('mmg.m : unavailable case')
end    
    
% Check configuration data and write original mesh in .msh format
if (nargin == 1)
    mshWriteMsh('original.msh',mesh);
elseif (nargin == 2)
    mshWriteMsh('original.msh',mesh,hnod);
elseif (nargin == 3)
    mshWriteMsh('original.msh',mesh);
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
file = '-in original.msh -out refined.msh ';  % mesh temporary files 

% Convert input to command
if (nargin==1)
    opt = ['-hsiz ',num2str(hsiz),' '];  
elseif (nargin==2)
    opt = ' ';  
elseif (nargin==3)
    opt = ['-hmin  ',num2str(hmin),' -hmax  ',num2str(hmax),' '];   
else
    error('mmg.m : unavailable case');
end
    
% Define verbose level : 0 (no infos), 1 default to 6 (everything)
info = ['-v ',num2str(0),' ']; 

% Execute mmg binaries
command = [os,bin,file,opt,info];
system(command);  

% Read refined mesh
[vtx,elt,val] = mshReadMsh('refined.msh');

% Update mesh
vtx  = nrm*vtx + ones(size(vtx,1),1)*ctr;
val  = nrm*val;
mesh = msh(vtx,elt);

% Clean meshes
delete('original.msh')
delete('refined.msh')

% Back to current directory
cd(here)
end
