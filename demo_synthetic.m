clc; clear; close all;
rng('default');
addpath([pwd,'/methods/GS']);
addpath([pwd,'/methods/CNMF']);
addpath([pwd,'/data_generation']);
% ============================================
load('S.mat');
load('endmember.mat');
load('HS_spec.mat');
load('MS_spec.mat');
abun_Cuprite = reshape(S_curr',480,480,[]);
W1 = 120; W2 = 120; L = W1*W2;
patchNum = 8;
M = 224;
N = 30;
% --------generate data
maxIteration = 100;
stopcriterion = 1e-5;
dsRatio = 4;
GauSigma = 1.7;
kernelSize = 11;
[G,B,S_LRSR] = Construct_Toeplitz_G(W1,W2,kernelSize,GauSigma,dsRatio);
F = Construct_F(HS_spec,MS_spec);
% ----CNMF settings
CNMF_opt.convthresh = stopcriterion;
CNMF_opt.maxIteraNum = maxIteration;
CNMF_opt.F = F;
CNMF_opt.G = G;
% -------set recorders
algorithmNum = 10;
trial = 100;
SNR_case = 3;
runtime = zeros(algorithmNum,trial,SNR_case);
PSNR = zeros(algorithmNum,trial,SNR_case);
SAM = zeros(algorithmNum,trial,SNR_case);
RMSE = zeros(algorithmNum,trial,SNR_case);
UIQI = zeros(algorithmNum,trial,SNR_case);
ERGAS = zeros(algorithmNum,trial,SNR_case);
for i = 1:trial
    % --------generate data
    W1_pos = randi(size(abun_Cuprite,1)-W1-1,1);
    W2_pos = randi(size(abun_Cuprite,2)-W2-1,1);
    S = abun_Cuprite(W1_pos:W1_pos+W1-1,W2_pos:W2_pos+W2-1,:);
    blk_sz1 = randfixedsumint(1,patchNum,W1);
    blk_sz2 = randfixedsumint(1,patchNum,W2);
    Y = zeros(W1,W2,M);
    for blk_num_1 = 1:patchNum
        for blk_num_2 = 1:patchNum
            A_blk = [endmember{1}(:,randperm(16,1)),endmember{2}(:,randperm(16,1)),endmember{3}(:,randperm(16,1)),...
                endmember{4}(:,randperm(16,1)),endmember{5}(:,randperm(16,1))];
            Y(sum(blk_sz1(1:blk_num_1-1))+1:sum(blk_sz1(1:blk_num_1-1))+blk_sz1(blk_num_1),sum(blk_sz2(1:blk_num_2-1))+1:sum(blk_sz2(1:blk_num_2-1))+blk_sz2(blk_num_2),:) = ...
                reshape((A_blk*reshape(S(sum(blk_sz1(1:blk_num_1-1))+1:sum(blk_sz1(1:blk_num_1-1))+blk_sz1(blk_num_1),sum(blk_sz2(1:blk_num_2-1))+1:sum(blk_sz2(1:blk_num_2-1))+blk_sz2(blk_num_2),:),[],5)')',blk_sz1(blk_num_1),[],M);
        end
    end
    dataset = Y;
    Y = reshape(Y,L,M)';
    Y_M = F*Y; Y_H = Y*G;
    V_H = randn(size(Y_H)); V_M = randn(size(Y_M));
    % --------different noisy scenarios
    SNR = 25;
    % ----observation generation
    YM_sigma = sqrt((sum(Y_M(:).^2)/(L))/(10^(SNR/10)))/sqrt(size(F,1));
    YM_noise = Y_M+YM_sigma*V_M;
    YH_sigma = sqrt((sum(Y_H(:).^2)/(L/dsRatio^2))/(10^(SNR/10)))/sqrt(size(F,2));
    YH_noise = Y_H+YH_sigma*V_H;
    HSI_noise = reshape(YH_noise',W1/dsRatio,W2/dsRatio,[]);
    MSI_noise = reshape(YM_noise',W1,W2,[]);
    % ----initialization
    A_init = SPA(YH_noise,N);
    S_init = rand(N,L);
    S_init = bsxfun(@rdivide,S_init,sum(S_init));
    X_init = A_init*S_init;
    % ---- GSA
    GSA_code = 1;
    tic;
    Y_GSA = GSA_wrapper(HSI_noise,MSI_noise,dsRatio);
    runtime(GSA_code,i) = toc;
    [PSNR(GSA_code,i),RMSE(GSA_code,i),ERGAS(GSA_code,i),SAM(GSA_code,i),UIQI(GSA_code,i),psnr_GSA,sam_GSA,mse_GSA] = ...
        quality_assessment(reshape(Y',W1,W2,[]),Y_GSA,0,1/dsRatio);
    fprintf('%gdB, trial %g: (GSA) time: %gs, PSNR: %g, SAM: %g.\n',...
        SNR,i,runtime(GSA_code,i),PSNR(GSA_code,i),SAM(GSA_code,i));
    % ---- CNMF
    CNMF_code = 2;
    tic;
    A_init = SPA(YH_noise,N);
    S_init = rand(N,L);
    S_init = bsxfun(@rdivide,S_init,sum(S_init));
    CNMF_opt.A_init = A_init;
    CNMF_opt.S_init = S_init;
    [A_CNMF,S_CNMF] = CNMF(HSI_noise,MSI_noise,CNMF_opt);
    runtime(CNMF_code,i) = toc;
    Y_CNMF = A_CNMF*S_CNMF;
    [PSNR(CNMF_code,i),RMSE(CNMF_code,i),ERGAS(CNMF_code,i),SAM(CNMF_code,i),UIQI(CNMF_code,i),psnr_CNMF,sam_CNMF,mse_CNMF] = ...
        quality_assessment(reshape(Y',W1,W2,[]),reshape(Y_CNMF',W1,W2,[]),0,1/dsRatio);
    fprintf('              (CNMF) time: %gs, PSNR: %g, SAM: %g.\n',...
        runtime(CNMF_code,i),PSNR(CNMF_code,i),SAM(CNMF_code,i));
    % ---- GLORIA
    GLORIA_code = 3;
    tic;
    Y_GLORIA_full = GLORIA_simplified(YH_noise,YM_noise,F,G,W1,W2,'MU',30/SNR);
    runtime(GLORIA_code,i) = toc;
    [PSNR(GLORIA_code,i),~,~,SAM(GLORIA_code,i)] = ...
        quality_assessment(reshape(Y',W1,W2,[]),reshape(Y_GLORIA_full',W1,W2,[]),0,1/dsRatio);
    fprintf('        (GLORIA-full) time: %gs, PSNR: %g, SAM: %g.\n',...
        runtime(GLORIA_code,i),PSNR(GLORIA_code,i),SAM(GLORIA_code,i));
end