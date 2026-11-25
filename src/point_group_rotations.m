function R = point_group_rotations(cif)
%POINT_GROUP_ROTATIONS  Proper rotations of the crystal in Cartesian coordinates.
%
%   R = point_group_rotations(cif)
%
%   Takes the rotation parts W of the symmetry operations of the CIF
%   (translations dropped), keeps those with det(W) = +1, converts them to
%   the Cartesian frame of CALCULATE_RECIPROCAL_LATTICE (a along x, b in
%   the xy plane) and removes duplicates.
%
%   R : 3 x 3 x N; N = 24 for cubic m-3m, 12 for 6/mmm, 8 for 4/mmm,
%       6 for -3m, 4 for mmm, 2 for 2/m, 1 for -1.

[~, ~, A] = calculate_reciprocal_lattice(struct('a', cif.cell(1:3), 'u', cif.cell(4:6)), [0 0 0], [1 0 0 0]);
M = A';                                   % columns = basis vectors a, b, c
R = zeros(3, 3, 0);
for k = 1:size(cif.W, 3)
    W = cif.W(:,:,k);
    if abs(det(W) - 1) > 1e-6, continue; end
    Rc = M * W / M;                       % fractional -> Cartesian
    dup = false;
    for j = 1:size(R, 3)
        if norm(R(:,:,j) - Rc) < 1e-6, dup = true; break; end
    end
    if ~dup, R(:,:,end+1) = Rc; end %#ok<AGROW>
end
end
