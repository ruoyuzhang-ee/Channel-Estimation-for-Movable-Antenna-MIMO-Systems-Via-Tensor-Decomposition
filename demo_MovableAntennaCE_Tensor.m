% 20240403
% Please refer to "Channel Estimation for Movable-Antenna MIMO Systems Via Tensor Decomposition," IEEE Wireless Communications Letters, 2024
% Ruoyu Zhang
 
addpath('tensorlab');

%% 1. Initialization
% This is the true number of antennas in movable antennas systems
N_BS_RF = 4;
N_MS_RF = 4;

c_light = 3*1e8;
fc = 3*1e9; 
%Virtual grid
N_BS_x = 40;
N_BS_y = 40;
N_BS = N_BS_x * N_BS_y;
N_MS_x = 40; 
N_MS_y = 40; 
N_MS = N_MS_x * N_MS_y;
lambda = c_light/fc;

Lt = 3;  %Number Channel Path
Lr = 3;  %Number Channel Path
Lt_est = Lt; 
Lr_est = Lr; 


dd_BS_x = 1/5 * lambda; 
dd_BS_y = 1/5 * lambda; 
dd_MS_x = 1/5 * lambda;
dd_MS_y = 1/5 * lambda; 
ArrayAperture_BS_x = N_BS_x * dd_BS_x;
ArrayAperture_BS_y = N_BS_y * dd_BS_y;
ArrayAperture_MS_x = N_MS_x * dd_MS_x;
ArrayAperture_MS_y = N_MS_y * dd_MS_y;


Mx = 20; My = 20;
M = Mx*My;
Nx = 20; Ny = 20;
N = Nx*Ny;
%Measurement grid
dd_tt_x = 1/5 * lambda;
dd_tt_y = 1/5 * lambda;
dd_rr_x = 1/5 * lambda;
dd_rr_y = 1/5 * lambda;


Num_Trial = 1000;
SNR_dB_list = [0:5:25];

%% mse
mse_H_ALL_Tensor = zeros(Num_Trial, length(SNR_dB_list));

for SNR_ii = 1:length(SNR_dB_list)
    sigma2 = 1;
    SNR_dB = SNR_dB_list(SNR_ii);
    Power = sigma2 * 10^(SNR_dB/10);
    
    
    for Trial_ii = 1:Num_Trial
        %% Channel
        load('H_BSx40_BSy40_MSx40_MSy40_Lt3_Lr3.mat');
        
        %Steering vector
        a_BS_matrix_x = exp(1j*2*pi * dd_BS_x/lambda *(0:N_BS_x-1).' .* AoD_x_true );
        a_BS_matrix_y = exp(1j*2*pi * dd_BS_y/lambda *(0:N_BS_y-1).' .* AoD_y_true );
        a_BS_matrix = kr(a_BS_matrix_x, a_BS_matrix_y);
        a_MS_matrix_x = exp(1j*2*pi * dd_MS_x/lambda *(0:N_MS_x-1).' .* AoA_x_true);
        a_MS_matrix_y = exp(1j*2*pi * dd_MS_y/lambda *(0:N_MS_y-1).' .* AoA_y_true);
        a_MS_matrix = kr(a_MS_matrix_x, a_MS_matrix_y);

        H = conj(a_MS_matrix) * Sigma_gain * a_BS_matrix.';

        
        
        %% Measurement Model
        %True Measurement Transmit Grid 
        a_BS_Mx = exp(1j*2*pi * dd_tt_x/lambda *(0:Mx-1).' .* AoD_x_true );
        a_BS_My = exp(1j*2*pi * dd_tt_y/lambda *(0:My-1).' .* AoD_y_true );
        a_BS_M = kr(a_BS_Mx, a_BS_My); 
        G = a_BS_M.';
        AA_alpha_t = G'; 
        %True Measurement Receive Grid 
        a_MS_Nx = exp(1j*2*pi * dd_rr_x/lambda *(0:Nx-1).' .* AoA_x_true );
        a_MS_Ny = exp(1j*2*pi * dd_rr_y/lambda *(0:Ny-1).' .* AoA_y_true );
        a_MS_N = kr(a_MS_Nx, a_MS_Ny); 
        BB_alpha_r = conj(a_MS_N);

        %Receive signal model
        %BS move, MS receive, for AoD estimatin
        initalPosition_MS = [1  Ny  (Nx-1)*Ny+1  Nx*Ny]; 
        a_MS_N_inital = a_MS_N(initalPosition_MS,:);
        noise = 1/sqrt(2) * (random('norm',0,1, M, N_MS_RF) + 1i*random('norm',0,1, M, N_MS_RF)); 
        noise_t = sqrt(sigma2) * noise;
        yy_t = sqrt(Power) * AA_alpha_t * Sigma_gain' * a_MS_N_inital.' + noise_t; 
        %MS move, BS receive  or MS move, BS transmit
        initalPosition_BS = [1  My  (Mx-1)*My+1  Mx*My]; 
        a_BS_M_inital = a_BS_M(initalPosition_BS,:);
        noise = 1/sqrt(2) * (random('norm',0,1, N, N_BS_RF) + 1i*random('norm',0,1, N, N_BS_RF)); 
        noise_r = sqrt(sigma2) * noise;
        yy_r = sqrt(Power) * BB_alpha_r * Sigma_gain * a_BS_M_inital.' + noise_r;
        
        %Receive signal model Tensor Formulation
        yy_t_Tensor = zeros(My,Mx,N_MS_RF);
        for nn = 1:N_MS_RF
            yy_t_Tensor(:,:,nn) = reshape(yy_t(:,nn), My, Mx);
        end
        yy_r_Tensor = zeros(Ny,Nx,N_BS_RF);
        for nn = 1:N_BS_RF
            yy_r_Tensor(:,:,nn) = reshape(yy_r(:,nn), Ny, Nx);
        end
        
        SysPara.N_BS_x = N_BS_x;
        SysPara.N_BS_y = N_BS_y;
        SysPara.N_MS_x = N_MS_x;
        SysPara.N_MS_y = N_MS_y;
        SysPara.dd_BS_x = dd_BS_x;
        SysPara.dd_BS_y = dd_BS_y;
        SysPara.dd_MS_x = dd_MS_x;
        SysPara.dd_MS_y = dd_MS_y;
        SysPara.Mx = Mx;
        SysPara.My = My;
        SysPara.Nx = Nx;
        SysPara.Ny = Ny;
        SysPara.dd_tt_x = dd_tt_x;
        SysPara.dd_tt_y = dd_tt_y;
        SysPara.dd_rr_x = dd_rr_x;
        SysPara.dd_rr_y = dd_rr_y;
        SysPara.lambda = lambda;
        SysPara.Power = Power;
        SysPara.initalPosition_MS = initalPosition_MS;
        SysPara.initalPosition_BS = initalPosition_BS;
        SysPara.Lr_est = Lr_est;
        SysPara.Lt_est = Lt_est;
        SysPara.yy_t = yy_t;
        SysPara.yy_r = yy_r;
        SysPara.yy_t_Tensor = yy_t_Tensor;
        SysPara.yy_r_Tensor = yy_r_Tensor;
        

        %% Tensor Decomposition
        [U_est] = cpd(yy_t_Tensor,Lt_est);
        BBhat1 = U_est{1,1};
        BBhat2 = U_est{1,2};
        BBhat3 = U_est{1,3};
        CChat1_AoD = BBhat1(:,1:Lt_est);
        CChat2_AoD = BBhat2(:,1:Lt_est);
        CChat3_AoD = BBhat3(:,1:Lt_est);

        [U_est] = cpd(yy_r_Tensor,Lr_est);
        BBhat1 = U_est{1,1};
        BBhat2 = U_est{1,2};
        BBhat3 = U_est{1,3};
        CChat1_AoA = BBhat1(:,1:Lr_est);
        CChat2_AoA = BBhat2(:,1:Lr_est);
        CChat3_AoA = BBhat3(:,1:Lr_est);
        
        %Est angle
        [ AoD_y_est ] = func_1D_estAngle_ESPRIT( CChat1_AoD, My, dd_tt_y, lambda);
        [ AoD_x_est ] = func_1D_estAngle_ESPRIT( CChat2_AoD, Mx, dd_tt_x, lambda);
        [ AoA_y_est ] = func_1D_estAngle_ESPRIT( CChat1_AoA, Ny, dd_rr_y, lambda);
        [ AoA_x_est ] = func_1D_estAngle_ESPRIT( CChat2_AoA, Nx, dd_rr_x, lambda);
        
        %Est channel
        EstPara_Tensor.AoD_y_est = AoD_y_est;
        EstPara_Tensor.AoD_x_est = AoD_x_est;
        EstPara_Tensor.AoA_y_est = AoA_y_est;
        EstPara_Tensor.AoA_x_est = AoA_x_est;
        [ outPara_Tensor ] = func_ChannelMatrixEst( SysPara, EstPara_Tensor );
        H_est = outPara_Tensor.H_est;
        
        nmse_H = norm(H_est - H,'fro')^2 / norm(H,'fro')^2;
        mse_H_ALL_Tensor(Trial_ii, SNR_ii) = nmse_H;

        
    end

end
figure; 
semilogy(SNR_dB_list, mean(mse_H_ALL_Tensor,1),'-*');  hold on; 
xlabel('SNR (dB)'); ylabel('NMSE (H)'); 



