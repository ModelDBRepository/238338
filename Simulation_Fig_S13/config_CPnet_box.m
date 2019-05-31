% cfg = config_CPnet_box([xdim, ydim], radius, hard_or_smooth)
% prepares a box-shaped (xdim x ydim) CP network where all the cells within
% `radius` are mutually connected.
%
% Input:
%   [xdim, ydim]   : x and y dimensions
%   radius         : radius for making connections
%   hard_or_smooth : scheme for making connections
%
% Output:
%   cfg    : geometric configuration
%
% Written by Sungho Hong, Computational Neuroscience, OIST, Japan
% Correspondence: Sungho Hong (shhong@oist.jp)
%
% February 9, 2018
%
function cfg = config_CPnet_box(xydim, radius, scheme)

if nargin==2
    scheme = 'hard';
end

ncell = prod(xydim);
position = zeros(ncell, 2);
count = 1;
for y=0:(xydim(2)-1)
    for x=0:(xydim(1)-1)
        position(count, 1) = x;
        position(count, 2) = y;
        count = count+1;
    end
end

if strcmp(scheme, 'smooth')
    dd = pdist(position)-1;
    dw = exp(-dd.^2/(2*radius^2));
    dw(dd>=2*radius)=0;
else
    dw = pdist(position);
    dw = double(dw<=radius);
end


coupling = squareform(dw);
coupling = sparse(coupling);

cfg = vars2struct(ncell, ...
                  position,...
                  coupling);

end