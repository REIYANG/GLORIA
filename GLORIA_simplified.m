function X = GLORIA_simplified(Y_H,Y_M,F,G,W1,W2,varargin)
%=====================================================================
% Programmer: 
% Ruiyuan Wu, E-mail: reiyangimenm@gmail.com
% Date: Sep 26, 2019
% -------------------------------------------------------
% Reference: 
% R. Wu, W.-K. Ma, X. Fu, and Q. Li, ``Hyperspectral Super-Resolution via Global-Local
% Low-Rank Matrix Estimation", submitted to IEEE Transactions on Geoscience and Remote Sensing, 2019. 
%======================================================================
% An implementation of HiBCD (Hybrid Inexact BCD) scheme for plain CoSMF
% Y_HiBCD = GLORIA_simplified(Y_H,Y_M,F,G,W1,W2,varargin)
%======================================================================
%  Input
%  Y_H is the spectral-spatial matrix form of the HS image
%  Y_M is the spectral-spatial matrix form of the MS image
%  F is the spectral decimation matrix
%  G is the spatial decimation matrix
%  W1 is the height of the image
%  W2 is the width of the image
%  varargin (optional):
%   - MU: regularization parameter
%   - INITIALIZATION
%----------------------------------------------------------------------
%  Output
%  X is the M-L SR estimation
%========================================================================
maxIter = 100;
patch_num = 16;
p = 1/2;
X = rand(size(Y_H,1),size(Y_M,2));
tau = 1e0;
% --------read the optional parameters
if (nargin-length(varargin))~=6
    error('Wrong number of required parameters');
end
if rem(length(varargin),2)==1
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MU'
                mu = varargin{i+1};
            case 'INITIALIZATION'
                X = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
eps = tau;
W1_num = sqrt(patch_num);
W2_num = sqrt(patch_num);
window_s1 = W1/W1_num;
window_s2 = W2/W2_num;
P_sub = cell(patch_num,1);
P_index = cell(patch_num,1);
k = 1;
for i = 1:W1_num
    P_1 = sparse([zeros((i-1)*window_s1,window_s1);eye(window_s1);zeros((W1_num-i)*window_s1,window_s1)]);
    for j = 1:W2_num
        P_2 = sparse([zeros((j-1)*window_s2,window_s2);eye(window_s2);zeros((W2_num-j)*window_s2,window_s2)]);
        P_sub{k} = sparse(kron(P_1,P_2));
        P_index{k} = find(sum(P_sub{k},2)~=0);
        k = k+1;
    end
end
theta_G = eigs(G'*G,1);
FTF = F'*F;
% theta_F = eigs(FTF,1);
X_prev = X;
gamma_prev = 1;
% ----start iteration
for count = 1:maxIter
    gamma = (1+sqrt(1+4*gamma_prev^2))/2;
    H = (gamma+gamma_prev-1)/gamma*X-(gamma_prev-1)/gamma*X_prev;
    X_prev = X;
    gamma_prev = gamma;
    % ----global patch
    [V,D] = svd(X*X');
    D_Z = p*(diag(D)+eps).^((p-2)/2);
    Z_inv = bsxfun(@times,V,D_Z')*V';
    mu_gr_lips = norm(FTF+mu*Z_inv,2);
    mu_gr_grad = mu*Z_inv*H;
    % ----local patch
    lr_grad = zeros(size(X));
    lips_lr = zeros(patch_num,1);
    for P_count = 1:patch_num
        P_i = P_sub{P_count};
        X_i = X*P_i;
        H_i = H*P_i;
        [V_i,D_i] = svd(X_i*X_i');
        DZ_i = p*(diag(D_i)+eps).^((p-2)/2);
        Zinv_i = bsxfun(@times,V_i,DZ_i')*V_i';
        lips_lr(P_count) = max(DZ_i);
        lr_grad(:,P_index{P_count}) = Zinv_i*H_i;
    end
    mu_lr_grad = mu*lr_grad;
    mu_lr_lips = mu*max(lips_lr);
    lips = theta_G+mu_gr_lips+mu_lr_lips;
    nabla_H = (H*G-Y_H)*G'+F'*(F*H-Y_M)+mu_gr_grad+mu_lr_grad;
    X = min(1,max(H-1/lips*nabla_H,0));
end
end
