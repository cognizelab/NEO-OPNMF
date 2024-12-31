function [xr, yr] = get_rotation(x, y, angle, direction)
% 2D coordinate rotation
%   [xr, yr] = get_rotation(x, y, angle, direction) returns rotated coordinates
%   Input parameters:
%       x, y: Original coordinates
%       angle: Rotation angle (in degrees)
%       direction: Rotation direction ('ccw' counterclockwise[default], 'cw' clockwise)
%   Output parameters:
%       xr, yr: Rotated coordinates

% Note:
% When the coordinate axes are rotated by an angle θ, 
% the y-component of the new coordinates of a data point can be interpreted as 
% the projection length of the original vector in a specific direction.
% 
% [xr,yr] = get_rotation(x, y, 30, 'cw');
% [pl] = get_projection(x, y, 60);
% "yr == pl"

% Check number of input arguments
if nargin < 3
    error('At least x, y coordinates and rotation angle are required');
end

% Set default rotation direction to counterclockwise
if nargin < 4
    direction = 'ccw';
end

% Adjust angle sign based on rotation direction
if strcmpi(direction, 'cw')
    angle = -angle;
elseif ~strcmpi(direction, 'ccw')
    error('Invalid rotation direction, use "ccw" or "cw"');
end

% Convert angle from degrees to radians
theta = angle * pi / 180;

% Create rotation matrix
R = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];

% Ensure x and y have the same dimensions
if ~isequal(size(x), size(y))
    error('x and y must have the same dimensions');
end

% Store original shape
original_size = size(x);

% Convert coordinates to column vectors
x = x(:);
y = y(:);

% Perform rotation transformation
rotated_coords = [x y] * R;

% Extract rotated coordinates and restore original shape
xr = reshape(rotated_coords(:,1), original_size);
yr = reshape(rotated_coords(:,2), original_size);
end