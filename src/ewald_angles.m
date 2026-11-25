function [tth, beta1, beta2, phi1, phi2, missed] = ewald_angles(g, wavelength)
%EWALD_ANGLES  Where reciprocal vectors diffract when the sample rotates.
%
%   [tth, beta1, beta2, phi1, phi2, missed] = ewald_angles(g, wavelength)
%
%   Geometry: incident beam along +x, wave vector K = 2*pi/wavelength;
%   the sample rotates about the z axis by the angle phi. A vector g
%   diffracts when |k0 + R_z(phi) g| = K, i.e. at two rotation angles.
%
%   g          : N x 3 reciprocal vectors at phi = 0 (1/A, including 2*pi)
%   wavelength : (A)
%
%   tth          : 2theta (deg)
%   beta1, beta2 : azimuth on the detector of the two solutions (deg)
%   phi1, phi2   : rotation angle of the two solutions (deg)
%   missed       : true for vectors that never reach the Ewald sphere
%                  (removed from all other outputs)
%   (identical to Rotace.m of the original scripts)

K = 2 * pi / wavelength;
gx = g(:,1);  gy = g(:,2);  gz = g(:,3);
g2 = gx.^2 + gy.^2 + gz.^2;

radicand = 4 * K^2 * (gx.^2 + gy.^2) - g2.^2;   % >= 0 if g can reach the sphere
missed = radicand < 0;
gx = gx(~missed);  gy = gy(~missed);  gz = gz(~missed);  g2 = g2(~missed);

x  = -g2 / (2 * K);                              % in-plane position at diffraction
y1 = -sqrt(radicand(~missed)) / (2 * K);
y2 = -y1;

phi_g = mod(atan2d(gy, gx) + 360, 360);
phi1  = mod(mod(atan2d(y2, x) + 360, 360) - phi_g + 360, 360);
phi2  = mod(mod(atan2d(y1, x) + 360, 360) - phi_g + 360, 360);

tth   = acos(1 + x / K) * 180 / pi;
beta1 = mod(270 - atan2d(gz, -y2), 360) - 90;
beta2 = mod(270 - atan2d(gz, -y1), 360) - 90;
end
