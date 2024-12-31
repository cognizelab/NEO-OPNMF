function [pl] = get_projection(x, y, angle, direction)
%  Calculates vector projection length
%   [pl] = get_projection(x, y, angle, direction) returns projection length
%   Input parameters:
%       x, y: Original coordinates
%       angle: Projection angle (in degrees)
%       direction: Angle direction ('ccw' counterclockwise[default], 'cw' clockwise)
%   Output parameters:
%       pl: Projection length

% Check number of input arguments
if nargin < 3
    error('At least x, y coordinates and projection angle are required');
end

% Set default direction to counterclockwise
if nargin < 4
    direction = 'ccw';
end

% Adjust angle sign based on direction
if strcmpi(direction, 'cw')
    angle = -angle;
elseif ~strcmpi(direction, 'ccw')
    error('Invalid direction, use "ccw" or "cw"');
end

% Convert angle to radians
theta = angle * pi / 180;

% Create unit vector in the projection direction
unit_vector = [cos(theta), sin(theta)];

% Ensure x and y have the same dimensions
if ~isequal(size(x), size(y))
    error('x and y must have the same dimensions');
end

% Save original size
original_size = size(x);

% Convert to column vectors
x = x(:);
y = y(:);

% Calculate dot product with unit vector
pl = x * unit_vector(1) + y * unit_vector(2);

% Restore original shape
pl = reshape(pl, original_size);
end