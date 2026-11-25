function q = orientation_grid(mesh_degree, region)
%ORIENTATION_GRID  Trial orientations from a subdivided icosahedron.
%
%   q = orientation_grid(mesh_degree)            directions in the octant x, y, z >= 0
%   q = orientation_grid(mesh_degree, 'octant')  the same
%   q = orientation_grid(mesh_degree, R)         directions in the fundamental
%                                                region of the rotations R (3x3xN)
%
%   An icosahedron is subdivided mesh_degree times, which gives an almost
%   uniform set of directions n on the unit sphere (10*4^d + 2 points).
%   Each direction is converted to the quaternion of the shortest rotation
%   that brings n to [0 0 1]; the remaining rotation about [0 0 1] (the
%   sample rotation axis) is scanned separately by phi0.
%
%   Directions related by a proper rotation of the crystal give the same
%   diffraction pattern, so only one region of the sphere needs to be
%   tested. With R given, a direction n is kept if c.n >= c.(R_k n) for all
%   rotations R_k, where c is a fixed generic direction; this selects the
%   fundamental region around c for any point group (1/24 of the sphere for
%   cubic, 1/12 hexagonal, 1/4 orthorhombic, the whole sphere for triclinic).
%   The octant is sufficient for cubic, hexagonal and tetragonal crystals
%   in their standard settings only.

if nargin < 2, region = 'octant'; end
mesh = SubdivideSphericalMesh(IcosahedronMesh, mesh_degree);
[~, n] = GetMeshData(mesh);
if ischar(region) || isstring(region)
    n = n(all(n >= 0, 2), :);
else
    c = [0.12; 0.31; 0.94];  c = c / norm(c);
    keep = true(size(n, 1), 1);
    cn = n * c;
    for k = 1:size(region, 3)
        keep = keep & (cn >= (n * region(:,:,k)') * c - 1e-9);
    end
    n = n(keep, :);
end
q = short_quat(n, repmat([0 0 1], size(n, 1), 1));
end
