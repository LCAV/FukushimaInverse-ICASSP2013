
% This code and all associated files are the supplementary material to the paper
% M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
% The Fukushima inverse problem, ICASSP 2013
%
% 2013 (c) M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
% All the code is published under a CC-BY-SA 3.0 License
% For details about the license, refer to http://creativecommons.org/licenses/by-sa/3.0/
%   * For attribution of non-commercial reuse of this work, a similar notice to this one is sufficient
%   * For attribution of commercial reuse of this work, please contact us.
% 
% Contact: marta.martinez-camara@epfl.ch

%%%%%%%%%%%%%%%%%%%
% MANUAL SETTINGS %
%%%%%%%%%%%%%%%%%%%

%% add the path to CVX if it is required
%addpath('../cvx/');


%%%%%%%%%%%%%%%
% SETTING ENV %
%%%%%%%%%%%%%%%

clear; close all; clc, % tabula rasa
disp('Seting up the environment...'); % let the user know we are starting to set up the simulation

% add necessary path
addpath('./Data');
addpath('./Code');
addpath('./Code/opt_routines');

% setup cvx thingy
cvx_setup;  % set up cvx

% setting random number generator
rng('shuffle'); % set the random seed to the clock (independent results each time this script runs)


%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%

%% target condition number for matrix cleaning
cleanConditionNumber = 1305;

%% lambda used in L1 regularization when recovering source from real data
lambda = 1e-5;

% This will normalize the solutions that are in Becquerels per 3 hours slot
% into Giga-Becquerels per second
norm_to_GBqs = 3*60*60 * 1e9;

%% we only keep the columns corresponding to the
%% 5 first days after the accidents (1 col == 3 hours)
%% Thus 5 days x (8x3 hours) x 3 heights = 120 columns
cols = 1:120;

% The transport matrix has entries too small for CVX
% we scale the whole system by a large number
scaling = 10^16;

% main parameters of the simulation with synthetic measurements
N = 100;                                % number of iterations for averaging
lmbds = logspace(-40,2,30);             % range of lambdas to try
epslns = logspace(-40,2,30);            % weights for the second derivative term
snr = 0:10:50;                          % signal-to-noise ratios we want to examine, in dB

% parameters for sensitivity to a priori guess
% These were picked because they give the best reconstruction
% when the a priori solution is wrong
lambda_sens  = 0.12;
epsilon_sens = 3.56;


%%%%%%%%%%%%%%%%
% LOADING DATA %
%%%%%%%%%%%%%%%%

%% Load and clean the matrix
disp('Load matrix...');

load('matrixGFSXe.mat', 'matrix'); 
load('measXe.mat', 'measurements');
load('aPrioriSource.mat', 'XaTotalInt');
aPrioriSource = norm_to_GBqs*XaTotalInt(:);


%%%%%%%%%%%%%%%%%%%%%%%%
% PRODUCE CLEAN MATRIX %
%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean the matrix
disp('Clean matrix...');
[M, V] = matrixCleaning(matrix, measurements, cleanConditionNumber);

% plot figure 2. (without the map)

figure(2);
matrix_dB = 10*log10(max(abs(matrix), 1e-30));
M_dB = 10*log10(max(abs(M(:,cols)), 1e-30));
clim = [ min(matrix_dB(:)) max(matrix_dB(:)) ];

subplot(1,2,1);
imagesc(matrix_dB, clim);
title('Fig. 2A : Original transport matrix');
set(gca,'DataAspectRatio',[1 1 1])

subplot(1,2,2);
imagesc(M_dB, clim);
title('Fig. 2B : Cleaned transport matrix');
set(gca,'DataAspectRatio',[1 1 1])
colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%
% REAL SOURCE RECOVERY %
%%%%%%%%%%%%%%%%%%%%%%%%

%% Reconstruct from real data
disp('Reconstruct emmisions from real data...');
realSource = scaling*reconstructSourceL1Pos(scaling*M(:,cols), V, lambda);
    
% plot figure 4. (without the fancy illustrator thingy)
figure(4);

% the time slots are 3 hours
t = 0:3:(3*length(cols)/3-1);

subplot(3,1,1);
plot(t, realSource(3:3:end)/norm_to_GBqs, 'b');
title('Fig. 4 : Reconstruction of emission rates in [GBq/s] using the proposed algorithm.');
ylabel('300m-1000m height');

subplot(3,1,2);
plot(t, realSource(2:3:end)/norm_to_GBqs, 'y');
ylabel('50m-300m height');

subplot(3,1,3);
plot(t, realSource(1:3:end)/norm_to_GBqs, 'r');
ylabel('0m-50m height');
xlabel('time in slices of 3 hours');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENSITIVITY TO INITIAL CONDITIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct an incorrect a priori solution
disp('Sensitivity to initial guess in smoothed approach...');
A_wrong_ap = aPrioriSource(cols);
A_wrong_ap(18:22) = A_wrong_ap(61:65);%wrong a priori solution
A_wrong_ap(61:65) = 0;

%% Create a 10dB noisy synthetic measurement
bn = awgn(scaling*M(:,cols)*aPrioriSource(cols), 10,'measured'); % noisy measurements

% compute the estimate once with correct apriori solution
[x_right_ap, mse_dummy, l_dummy, e_dummy] = stohl_et_al_ap(scaling*M(:,cols), bn, aPrioriSource(cols), lambda_sens, epsilon_sens, aPrioriSource(cols));

% and once with wrong a priori solution
[x_wrong_ap, mse_dummy, l_dummy, e_dummy] = stohl_et_al_ap(scaling*M(:,cols), bn, aPrioriSource(cols), lambda_sens, epsilon_sens, A_wrong_ap);

% plot figure 1. from paper

figure(1);

t = 0:3:(120*3-1);
lim = [-1e5 3.5e5];

subplot(2,1,1);
plot(t, x_right_ap/norm_to_GBqs, 'r', t, aPrioriSource(cols)/norm_to_GBqs, 'g');
legend('Solution', 'A-priori guess');
ylabel('Using expert knowledge a-priori guess.');
%ylim(lim);
title('Solutions with Different A-priori Guesses');

subplot(2,1,2);
plot(t, x_wrong_ap/norm_to_GBqs, 'r', t, A_wrong_ap/norm_to_GBqs, 'g');
legend('Solution', 'A-priori guess');
ylabel('Using a-priori guess with modified main peak.');
%ylim(lim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHETIC MEASUREMENTS EXPERIMENT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run the simulation
disp('Synthetic measurements experiment...');
[MSE_L1, MSE_Tik, MSE_Der] = reconstructFromSyntheticData_Xe(M(:,cols), aPrioriSource(cols), N, lmbds, epslns, snr, scaling);

% plot figure 3.

figure(3);

plot(snr, MSE_L1, 'b.-', snr, MSE_Tik, 'r', snr, MSE_Der, 'g--');
legend('L1','Tikhonov','Smoothed');
title('Fig. 3 : SNR of recovery from synthetic measurements.');


