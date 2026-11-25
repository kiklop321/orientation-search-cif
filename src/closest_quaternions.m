function q_near = closest_quaternions(q_best, mesh_degree, region)
%CLOSEST_QUATERNIONS  Orientations of a finer grid sorted by distance to q_best.
%
%   q_near = closest_quaternions(q_best, mesh_degree, region)
%
%   Builds ORIENTATION_GRID(mesh_degree, region) and returns its quaternions
%   sorted by the squared Euclidean distance to q_best (closest first).

if nargin < 3, region = 'octant'; end
q = orientation_grid(mesh_degree, region);
[~, order] = sort(sum((q - q_best).^2, 2), 'ascend');
q_near = q(order, :);
end
