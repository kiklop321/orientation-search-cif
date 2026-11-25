function euler = euler_bunge(g)
%EULER_BUNGE  Bunge Euler angles of an orientation matrix.
%
%   euler = euler_bunge(g)
%
%   g     : 3x3 orientation matrix, sample -> crystal (rows = the crystal
%           axes a, b, c written in the sample frame); for a crystal -> sample
%           matrix L use g = L'
%   euler : [phi1 Phi phi2] (deg), Bunge (z-x'-z''), g = Rz(phi2) Rx(Phi) Rz(phi1)
%           (passive rotations):
%             Phi  = acos(g33)
%             phi1 = atan2(g31, -g32),  phi2 = atan2(g13, g23)
%           for Phi = 0 only phi1 + phi2 is defined; phi2 = 0 is then used

Phi = acosd(max(-1, min(1, g(3,3))));
if sind(Phi) > 1e-9
    phi1 = atan2d(g(3,1), -g(3,2));
    phi2 = atan2d(g(1,3), g(2,3));
else
    phi1 = atan2d(g(1,2), g(1,1));
    phi2 = 0;
end
euler = [mod(phi1, 360), Phi, mod(phi2, 360)];
euler(abs(euler) < 1e-9 | abs(euler - 360) < 1e-9) = 0;      % no -0, 3e-14 or 360
end
