function q = short_quat(v1, v2)
%SHORT_QUAT  Quaternions of the shortest rotations taking v1 to v2.
%
%   q = short_quat(v1, v2)
%
%   v1, v2 : N x 3 vectors (need not be normalised)
%   q      : N x 4 unit quaternions [w x y z]; row i rotates about the axis
%            v1(i,:) x v2(i,:) by the angle between the two vectors, so
%            that quatrotate(quatconj(q(i,:)), v1(i,:)) is parallel to v2(i,:).

v1n = v1 ./ vecnorm(v1, 2, 2);          % (not NORMALIZE: MTEX has its own on the path)
v2n = v2 ./ vecnorm(v2, 2, 2);
cos_th = dot(v1n, v2n, 2);
ax = cross(v1n, v2n, 2);
% parallel or antiparallel vectors: the cross product is zero; any axis
% perpendicular to v1 will do (for parallel ones the angle is 0 anyway)
for i = find(vecnorm(ax, 2, 2) < 1e-12)'
    ax(i, :) = cross(v1n(i, :), [1 0 0]);
    if norm(ax(i, :)) < 1e-6, ax(i, :) = cross(v1n(i, :), [0 1 0]); end   % v1 along x
end
axis_n = ax ./ vecnorm(ax, 2, 2);
cos_th = max(-1, min(1, cos_th));

q = [sqrt((1 + cos_th) / 2), sqrt((1 - cos_th) / 2) .* axis_n];
end
