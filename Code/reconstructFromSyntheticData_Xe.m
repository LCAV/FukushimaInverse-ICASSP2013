function [MSE_L1, MSE_Tik, MSE_Der] = reconstructFromSyntheticData_Xe(A, xArt, N, lmbds, epslns, snr, scaling)
% [MSE_L1, MSE_Tik, MSE_Der] = reconstructFromSyntheticData_Xe(A, xArt, N, lmbds, epslns, snr)
% 
% This functions runs simulation of emission recovery for synthetic measurements using different techniques
% and with additive white gaussian noise assumption.
% We compare three different methods:
%
% - L1 sparse regularization
% - Tikhonov regularization
% - Tikhonov with 2nd derivative smoothing
%
% Arguments
% A : transport matrix
% xArt : Synthetic emmision vector
% N : number of repetition of the simulations
% lmbds : a range of lambdas (regularization parameter) to test
% epslns : a range of epsilons (regularization parameter) to test
% snr : a range of SNR to try
% scaling : pre-conditioning for L1 regularization 
%
% Return
% MSE_L1 : Mean Squared Error of L1 sparse regularization
% MSE_Tik : Mean Squared Error of Tikhonov regularization
% MSE_Der : Mean Squared Error of Tikhonov with 2nd derivative smoothing


% This code and all associated files are the supplementary material to the paper
% M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
% The Fukushima inverse problem, ICASSP 2013
%
% 2013 (c) M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
% 
% All the code is published under a CC-BY-SA 3.0 License
% For details about the license, refer to http://creativecommons.org/licenses/by-sa/3.0/
%   * For attribution of non-commercial reuse of this work, a similar notice to this one is sufficient
%   * For attribution of commercial reuse of this work, please contact us.
% 
% Contact: marta.martinez-camara@epfl.ch


%% Make xArt column vector
xArt = xArt(:);
n = length(xArt);

%% empty containers
MSE_L1  = zeros(N, length(snr));
MSE_Tik = zeros(N, length(snr));
MSE_Der = zeros(N, length(snr));

Lambda_L1  = zeros(N, length(snr));
Lambda_Tik = zeros(N, length(snr));
Lambda_Der = zeros(N, length(snr));

Epsilon_Der = zeros(N, length(snr));

x_der = zeros(n,N, length(snr));
x_l1 = zeros(n,N, length(snr));
x_tik = zeros(n,N, length(snr));

% construct artificial measurements without noise
b = A*xArt; 

%% the actual experiments
timings = zeros(N,1); % record of how long each iteration took

for inst = 1:N % loop for mse instances
  disp(['Iteration ',num2str(inst),'...']); % show the iteration number
  tic; % start timing

  for snrInd = 1:size(snr,2) % loop through all the noise levels

    %% Create noisy source vector
    bn = awgn(b, snr(snrInd),'measured'); % noisy measurements

    % === OPTIMIZATION 1 === (L1)
    [x_l1(:,inst,snrInd), MSE_L1(inst, snrInd), Lambda_L1(inst, snrInd)] = l1min(A, bn, xArt, lmbds, scaling);
    
    % === OPTIMIZATION 2 === (Tikhonov)
    [x_tik(:,inst,snrInd), MSE_Tik(inst, snrInd), Lambda_Tik(inst, snrInd)] = tikhonov(A, bn, xArt, lmbds);
    
    % === OPTIMIZATION 3 === (Tikhonov with 2nd derivative smoothing)
    [x_der(:,inst,snrInd), MSE_Der(inst, snrInd), Lambda_Der(inst, snrInd), Epsilon_Der(inst, snrInd)] = stohl_et_al(A, bn, xArt, lmbds, epslns);
  end % for loop snr

  timings(inst) = toc;
  disp(['Iteration ', num2str(inst), ' : ', num2str(timings(inst)), ' sec']);
end % for loop mse

% average the error
normal = mean(xArt.^2);
MSE_L1 = 10*log10(normal./mean(MSE_L1,1));
MSE_Tik = 10*log10(normal./mean(MSE_Tik,1));
MSE_Der = 10*log10(normal./mean(MSE_Der,1));

