function [pos, elements] = expand_positions(cif)
%EXPAND_POSITIONS  All atoms of the unit cell from the asymmetric unit of a CIF.
%
%   [pos, elements] = expand_positions(cif)
%
%   cif      : struct from READ_CIF
%   pos      : M x 5 [x y z occupancy type] with x, y, z in [0, 1);
%              type indexes elements
%   elements : cell array of the different elements
%
%   Every atom of cif.atoms is transformed by every symmetry operation,
%   mapped into the unit cell, and duplicates (positions closer than 1e-3
%   in fractional coordinates, taking the periodicity into account) are
%   removed.

[elements, ~, type] = unique(cif.atoms.element, 'stable');
pos = zeros(0, 5);
tol = 1e-3;
% Coordinates written with few digits (0.3333, 0.6667, 0.1667, ...) are
% put on the nearest multiple of 1/24 if they are closer than 5e-4 to it,
% so that symmetry-related atoms cancel exactly in forbidden reflections.
snap = @(x) (abs(x - round(x * 24) / 24) < 5e-4) .* (round(x * 24) / 24) + ...
            (abs(x - round(x * 24) / 24) >= 5e-4) .* x;
for a = 1:height(cif.atoms)
    r = snap([cif.atoms.x(a); cif.atoms.y(a); cif.atoms.z(a)]);
    P = zeros(size(cif.W, 3), 3);
    for k = 1:size(cif.W, 3)
        P(k, :) = (cif.W(:,:,k) * r + cif.w(:,k))';
    end
    P = mod(P, 1);
    P(abs(P - 1) < tol) = 0;
    keep = true(size(P, 1), 1);
    for k = 2:size(P, 1)                       % remove symmetry-equivalent duplicates
        dd = abs(P(1:k-1, :) - P(k, :));
        dd = min(dd, 1 - dd);
        if any(all(dd < tol, 2) & keep(1:k-1)), keep(k) = false; end
    end
    P = P(keep, :);
    pos = [pos; P, repmat([cif.atoms.occupancy(a), type(a)], size(P, 1), 1)]; %#ok<AGROW>
end
end
