function [ AoD_y_est ] = func_1D_estAngle_ESPRIT( CChat1, My, dd_tt_y, lambda)

[~,L] = size(CChat1);
AoD_y_est = zeros(1,L);

for iiL = 1: L
    aa = CChat1(1:My-1,iiL);
    bb = CChat1(2:My,iiL);
    zz_est = pinv(aa) * bb;
    phase_est = atan2(imag(zz_est),real(zz_est));
    AoD_y_est(iiL) = - phase_est/ (2*pi * dd_tt_y) * lambda;

    if AoD_y_est(iiL) < -1
        AoD_y_est(iiL) = -1;  %AoD_y_est(iiL) = max(angle, -1);
    elseif AoD_y_est(iiL) > 1
        AoD_y_est(iiL) = 1;  %AoD_y_est(iiL) = min(angle, 1);
    end
end

    
end