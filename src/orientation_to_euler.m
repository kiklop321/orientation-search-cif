function [euler, g] = orientation_to_euler(q, phi0)
%ORIENTATION_TO_EULER  Bunge Euler angles of an orientation found by the search.
%
%   [euler, g] = orientation_to_euler(q, phi0)
%
%   q, phi0 : orientation from SEARCH_ORIENTATION / orientation.mat
%   euler   : [phi1 Phi phi2] (deg), Bunge convention (z-x'-z''), as used in
%             EBSD (e.g. MTEX, OIM, AZtec)
%   g       : 3x3 orientation matrix, sample -> crystal (rows = the crystal
%             axes a, b, c written in the sample frame)
%
%   Sample frame = the laboratory in the middle of the exposure:
%   X = beam direction, Z = rotation axis of the sample, Y = Z x X.
%   (EBSD software usually uses another sample frame, e.g. Z = normal of the
%   sample surface; then the frames must first be related.)
%
%   Steps:
%   1. quaternion -> rotation matrix:  R = quat2dcm(q)  (= q as a 3x3 matrix)
%   2. add phi0 (a rotation of the sample about Z):
%          L = Rz(phi0) * R'      crystal -> sample (columns = a, b, c in the sample)
%          g = L'                 sample -> crystal (Bunge's orientation matrix)
%   3. g = Rz(phi2) * Rx(Phi) * Rz(phi1) (passive rotations), i.e.
%          g(3,3) = cos(Phi)
%          g(3,1) = sin(phi1) sin(Phi),   g(3,2) = -cos(phi1) sin(Phi)
%          g(1,3) = sin(phi2) sin(Phi),   g(2,3) =  cos(phi2) sin(Phi)
%      so  Phi = acos(g33),  phi1 = atan2(g31, -g32),  phi2 = atan2(g13, g23)
%      (for Phi = 0 only phi1 + phi2 is defined; phi2 = 0 is then used).
%
%   Because q only tilts the crystal (q(4) = 0) and phi0 only turns it about
%   Z, Phi and phi2 come from q alone, and phi0 simply adds to phi1.
%   The crystal symmetry gives equivalent angle sets (24 for cubic); this is
%   the one belonging to the q found by the search.

Rz = [cosd(phi0) -sind(phi0) 0; sind(phi0) cosd(phi0) 0; 0 0 1];
L = Rz * quat2dcm(q(:)' / norm(q))';
g = L';

euler = euler_bunge(g);                 % step 3, see EULER_BUNGE
end
