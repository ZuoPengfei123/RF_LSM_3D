function varargout = obtain_crust1_QB(varargin)
lat = varargin{1};
lon = varargin{2};
H = varargin{3};
if_topo = varargin{4};
if nargin < 3
    H = 0;
end
if nargin < 4
    if_topo = 0;
end
nsedi=0;
%------------------2024/8/18------------------------------------

model = [2.8   2.7  1.5  1.4200;
    8.0  5.29000    3.2100    2.4600; 
    18.0  6.1000    3.6000    2.8100;
    23.0  6.48000    3.7400    2.9200;
    40.0  6.63       3.79      2.9500;
    45.0   6.72       3.80      2.9800;
    59.0  8.200    4.4300    3.3600];

varargout{1} = model;
varargout{2} = nsedi;
end