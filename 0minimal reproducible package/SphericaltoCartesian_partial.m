function [dr_dx, dr_dy, dr_dz, dtheta_dx, dtheta_dy, dtheta_dz, dphi_dx, dphi_dy, dphi_dz] ...
                    = SphericaltoCartesian_partial(r,theta,phi)
%%
% Obtain the partial derivatives of spherical coordinates with respect to
% Cartesian coordinates.
% r--radial coordinate, should be a single number
% theta--polar angle, radian mesurement (SHOULD be row vector)   [0,pi]
% phi--azimuthal angle, radian mesurement (SHOULD be row vector) [0,2*pi]
%%

if r == 0
    error('r should not be zero!');
end
if theta == 0
    theta = theta + eps;
end

dr_dx = cos(phi') * sin(theta);
dr_dy = sin(phi') * sin(theta);
dr_dz = ones(length(phi), 1) * cos(theta);
dtheta_dx = cos(phi') * cos(theta) / r;
dtheta_dy = sin(phi') * cos(theta) / r;
dtheta_dz = - ones(length(phi), 1) * sin(theta) / r;
dphi_dx = - sin(phi') * sin(theta).^(-1) / r;
dphi_dy = cos(phi') * sin(theta).^(-1) / r;
dphi_dz = 0;
