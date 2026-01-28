function [ outPara ] = func_ChannelMatrixEst( SysPara, EstPara )

N_BS_x = SysPara.N_BS_x;
N_BS_y = SysPara.N_BS_y;
N_MS_x = SysPara.N_MS_x;
N_MS_y = SysPara.N_MS_y;
dd_BS_x = SysPara.dd_BS_x;
dd_BS_y = SysPara.dd_BS_y;
dd_MS_x = SysPara.dd_MS_x;
dd_MS_y = SysPara.dd_MS_y;
Mx = SysPara.Mx;
My = SysPara.My;
Nx = SysPara.Nx;
Ny = SysPara.Ny;
dd_tt_x = SysPara.dd_tt_x;
dd_tt_y = SysPara.dd_tt_y;
dd_rr_x = SysPara.dd_rr_x;
dd_rr_y = SysPara.dd_rr_y;
lambda = SysPara.lambda;
Power = SysPara.Power;
initalPosition_MS = SysPara.initalPosition_MS;
initalPosition_BS = SysPara.initalPosition_BS;
Lr_est = SysPara.Lr_est;
Lt_est = SysPara.Lt_est;
yy_t = SysPara.yy_t;
yy_r = SysPara.yy_r;
yy_t_Tensor = SysPara.yy_t_Tensor;
yy_r_Tensor = SysPara.yy_r_Tensor;

AoD_y_est = EstPara.AoD_y_est;
AoD_x_est = EstPara.AoD_x_est;
AoA_y_est = EstPara.AoA_y_est;
AoA_x_est = EstPara.AoA_x_est;

%% Measurement matrix based on the est angles
%Model Transmit Est
a_BS_Mx_est = exp(1j*2*pi * dd_tt_x/lambda *(0:Mx-1).' .* AoD_x_est );
a_BS_My_est = exp(1j*2*pi * dd_tt_y/lambda *(0:My-1).' .* AoD_y_est );
a_BS_M_est = kr(a_BS_Mx_est, a_BS_My_est);  %等价(6)的 G^{T}
G_est = a_BS_M_est.';
AA_alpha_t_est = G_est'; %conj(a_BS_M);
%Model Receive Est
a_MS_Nx_est = exp(1j*2*pi * dd_rr_x/lambda *(0:Nx-1).' .* AoA_x_est );
a_MS_Ny_est = exp(1j*2*pi * dd_rr_y/lambda *(0:Ny-1).' .* AoA_y_est );
a_MS_N_est = kr(a_MS_Nx_est, a_MS_Ny_est);  %等价(9)的 B
BB_alpha_r_est = conj(a_MS_N_est);

%yy_t = sqrt(Power) * AA_alpha_t * Sigma_gain' * a_MS_N_inital.' + noise_t;
%yy_t' = sqrt(Power) * conj(a_MS_N_inital) * Sigma_gain' * AA_alpha_t' + noise_t;
Phi_t_est = sqrt(Power) * kron( conj(AA_alpha_t_est), conj(a_MS_N_est(initalPosition_MS,:)) );
Phi_r_est = sqrt(Power) * kron( a_BS_M_est(initalPosition_BS,:), BB_alpha_r_est );

yy = [vec(yy_t'); vec(yy_r)];
Phi = [Phi_t_est; Phi_r_est];
gamma_est = pinv(Phi) * yy;
Sigma_gain_est = reshape(gamma_est, Lr_est, Lt_est);

a_BS_matrix_x_est = exp(1j*2*pi * dd_BS_x/lambda *(0:N_BS_x-1).' .* AoD_x_est );
a_BS_matrix_y_est = exp(1j*2*pi * dd_BS_y/lambda *(0:N_BS_y-1).' .* AoD_y_est );
a_BS_matrix_est = kr(a_BS_matrix_x_est, a_BS_matrix_y_est);
a_MS_matrix_x_est = exp(1j*2*pi * dd_MS_x/lambda *(0:N_MS_x-1).' .* AoA_x_est);
a_MS_matrix_y_est = exp(1j*2*pi * dd_MS_y/lambda *(0:N_MS_y-1).' .* AoA_y_est);
a_MS_matrix_est = kr(a_MS_matrix_x_est, a_MS_matrix_y_est);

H_est = conj(a_MS_matrix_est) * Sigma_gain_est * a_BS_matrix_est.';

outPara.Sigma_gain_est = Sigma_gain_est;
outPara.H_est = H_est;
end

