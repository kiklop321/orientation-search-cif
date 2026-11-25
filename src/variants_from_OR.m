function V = variants_from_OR(parent, child, OR)
%VARIANTS_FROM_OR  All crystallographic variants of a product phase in a parent crystal.
%
%   V = variants_from_OR(parent, child, OR)
%
%   parent, child : structs from READ_CIF (parent: the crystal whose
%                   orientation is searched, e.g. the beta grain; child: the
%                   phase that formed in it, e.g. alpha'')
%   OR            : orientation relationship, two pairs of parallel directions
%                   OR.parent_dirs = [u1 v1 w1; u2 v2 w2]  (indices in the parent cell)
%                   OR.child_dirs  = [u1 v1 w1; u2 v2 w2]  (indices in the child cell)
%                   row 1 of parent_dirs || row 1 of child_dirs, and row 2 || row 2;
%                   the two directions of each pair must be perpendicular.
%
%   V : 3 x 3 x N rotation matrices; V(:,:,j) turns the Cartesian frame of
%       the child into the frame of the parent for variant j.
%
%   Variant 1 follows directly from the OR. Applying all proper rotations of
%   the parent gives the other variants; two of them are the same variant
%   when they differ only by a rotation of the child itself
%   (V_i = V_j * R_child). The number of variants is
%   |rotations of parent| / |rotations common to both|, e.g. 12 for
%   beta (m-3m) -> alpha'' (Cmcm) with the Burgers-type relationship.

frame_p = frame(parent.cell, OR.parent_dirs);
frame_c = frame(child.cell,  OR.child_dirs);
V1 = frame_p / frame_c;                         % child Cartesian -> parent Cartesian

Rp = point_group_rotations(parent);
Rc = point_group_rotations(child);
V = zeros(3, 3, 0);
for k = 1:size(Rp, 3)
    Vk = Rp(:,:,k) * V1;
    known = false;
    for j = 1:size(V, 3)
        for m = 1:size(Rc, 3)
            if norm(Vk - V(:,:,j) * Rc(:,:,m)) < 1e-6, known = true; break; end
        end
        if known, break; end
    end
    if ~known, V(:,:,end+1) = Vk; end %#ok<AGROW>
end
end

function F = frame(cell, dirs)
% orthonormal frame (columns) built from two perpendicular lattice directions
[~, ~, A] = calculate_reciprocal_lattice(struct('a', cell(1:3), 'u', cell(4:6)), [0 0 0], [1 0 0 0]);
d1 = A' * dirs(1,:)';  d1 = d1 / norm(d1);
d2 = A' * dirs(2,:)';  d2 = d2 / norm(d2);
if abs(dot(d1, d2)) > 1e-3
    error('variants_from_OR:notPerpendicular', ...
        'The two directions [%s] and [%s] are not perpendicular.', num2str(dirs(1,:)), num2str(dirs(2,:)));
end
F = [d1, d2, cross(d1, d2)];
end
