function [h, B, A, V] = calculate_reciprocal_lattice(lattice, hkl, rot)
%CALCULATE_RECIPROCAL_LATTICE  Cartesian reciprocal-lattice vectors.
%
%   [h, B, A, V] = calculate_reciprocal_lattice(lattice, hkl, rot)
%
%   Input
%     lattice.a : lattice parameters [a b c] (A)
%     lattice.u : lattice angles [alpha beta gamma] (deg)
%     hkl       : N x 3 Miller indices
%     rot       : orientation quaternion [w x y z] applied to the direct
%                 lattice; [1 0 0 0] = no rotation
%
%   Output
%     h : N x 3 reciprocal-lattice vectors, h = hkl * B (1/A, incl. 2*pi)
%     B : 3 x 3 reciprocal basis (rows b1*, b2*, b3*)
%     A : 3 x 3 direct basis (rows a1, a2, a3); a1 along x, a2 in xy plane
%     V : unit-cell volume (A^3)

u = lattice.u * pi / 180;
a = lattice.a;

% Direct basis in Cartesian coordinates
A = zeros(3);
A(1,:) = [a(1), 0, 0];
A(2,:) = [a(2) * cos(u(3)), a(2) * sin(u(3)), 0];
A(3,1) = a(3) * cos(u(2));
A(3,2) = a(3) * (cos(u(1)) - cos(u(2)) * cos(u(3))) / sin(u(3));
A(3,3) = sqrt(a(3)^2 - A(3,1)^2 - A(3,2)^2);
V = det(A);

A = quatrotate(quatconj(rot), A);

% Reciprocal basis
B = 2 * pi * inv(A)';
h = hkl * B;
end
