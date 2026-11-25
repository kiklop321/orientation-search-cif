function [g, hkl_all, variant, F2] = phase_reflections(phase, wavelength)
%PHASE_REFLECTIONS  Reciprocal-lattice vectors of the allowed reflections of a phase.
%
%   [g, hkl_all, variant, F2] = phase_reflections(phase, wavelength)
%
%   phase      : struct from PHASE_FROM_CIF
%   wavelength : X-ray wavelength (A)
%
%   All (h k l) with |h|,|k|,|l| <= phase.hkl_max ('auto' = all that can
%   diffract at this wavelength) and d >= wavelength/2
%   (reachable, 2theta <= 180 deg) are generated and their structure factor
%       F = sum_n occ_n Z_n exp(2 pi i (h x_n + k y_n + l z_n))
%   is computed with the atomic number Z as scattering factor. Reflections
%   with |F|^2 < 1e-8 |F(000)|^2 (systematic extinctions) and 000 are
%   removed. For deciding which reflections exist this is equivalent to the
%   full calculation with angle-dependent form factors, unless atoms of
%   different elements cancel accidentally.
%
%   g       : M x 3 reciprocal vectors (1/A, including 2*pi) in the parent
%             frame, for every variant
%   hkl_all : M x 3 Miller indices
%   variant : M x 1 variant number
%   F2      : M x 1 |F|^2 (electrons^2, f = Z) of each reflection

n = phase.hkl_max;
if ischar(n) || isstring(n) || isinf(n)
    % 'auto': every index that can reach the Ewald sphere (|g| <= 2K, i.e.
    % d >= wavelength/2); the largest one is about 2*max(a, b, c)/wavelength
    n = ceil(2 * max(phase.cell(1:3)) / wavelength);
end
[h, k, l] = meshgrid(-n:n);
hkl = [h(:) k(:) l(:)];
hkl = hkl(any(hkl ~= 0, 2), :);

lat = struct('a', phase.cell(1:3), 'u', phase.cell(4:6));
g0 = calculate_reciprocal_lattice(lat, hkl, [1 0 0 0]);
reach = vecnorm(g0, 2, 2) <= 2 * (2 * pi / wavelength);   % |g| <= 2K
hkl = hkl(reach, :);  g0 = g0(reach, :);

Zat = phase.Z(phase.positions(:,5));
w = phase.positions(:,4) .* Zat(:);                       % occupancy * Z of every atom
F = exp(2i * pi * hkl * phase.positions(:,1:3)') * w;
F000 = sum(w);
keep = abs(F).^2 >= 1e-8 * F000^2;
hkl = hkl(keep, :);  g0 = g0(keep, :);  F2one = abs(F(keep)).^2;

nv = size(phase.variants, 3);
nh = size(hkl, 1);
g = zeros(nv * nh, 3);  hkl_all = zeros(nv * nh, 3);  variant = zeros(nv * nh, 1);  F2 = zeros(nv * nh, 1);
for j = 1:nv
    rows = (j - 1) * nh + (1:nh);
    g(rows, :) = (phase.variants(:, :, j) * g0')';
    hkl_all(rows, :) = hkl;
    variant(rows) = j;
    F2(rows) = F2one;
end
end
