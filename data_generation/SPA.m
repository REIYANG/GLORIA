%===================================================================================================
% Use successive projection algorithm to identify endmembers
%---------------------------------------------------------------------------------------------------
% Reference: 
% W.-K. Ma, J. M. Bioucas-Dias, T.-H. Chan, N. Gillis, P. Gader, A. Plaza, A. Ambikapathi 
% and C.-Y. Chi, "A signal processing perspective on hyperspectral unmixing," IEEE Signal
% Processing Magazine, vol. 31, no. 1, pp. 67-81, Jan 2014. 
%---------------------------------------------------------------------------------------------------
%------------------------input-------------------------
% Y is the observation matrix.
% N is the number of endmember.
%------------------------output------------------------
% A_est is the estimated endmember.
%===================================================================================================

function [A_est] = SPA(Y,N)

M = size(Y,1); % the spectral band
P = eye(M); % initialize the projection with unit matrix
A_est = zeros(M,N); % initialize the endmember matrix

for i = 1:1:N
    [~,index] = max(sum((P*Y).*(P*Y)));
    A_est(:,i) = Y(:,index);
    P = P-(P*Y(:,index))*((P*Y(:,index))'*P)/norm(P*Y(:,index),2)^2; % update the orthogonal projector
end