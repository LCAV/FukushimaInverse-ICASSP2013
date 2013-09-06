function [xhat, err, lambdastar, epsilonstar] = stohl_et_al(A, bn, s, lambdas, epsilons)
% Solves the inverse problem using Tikhonov regularization
%     min |A*x - bn|_2^2 + lambda*|x|_2^2 + epsilon*|D2*x|_2^2
% for all given values of lambda.
%
% Param:
%   A       : forward matrix
%   bn      : noisy measurement vector
%   s       : original source vector (bn = A*s+noise)
%   lambdas : values of the optimization parameter
%
% Return the optimal solution, the (mean squared) error of that solution,
% and the optimal value of lambda.


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


% length of solution
n = length(s);

%% 2nd derivative matrix
row2=zeros(n,1); % 1st row of D2
row2(1)=1;
col2=zeros(n,1); % 1st col of D2
col2(1)=1;
col2(2)=-2;
col2(3)=1;
D2=toeplitz(col2,row2); % second derivative operator

for i=1:length(lambdas)
  for j=1:length(epsilons) % loop through all the epsilons
    %Analytic solution
    x=(A'*A+lambdas(i)*eye(n)+epsilons(j)*(D2'*D2))\(A'*bn);

    % Check for minimum over lambdas
    mse=mean((x-s).^2);
    if (i == 1 || mse < err)
      err = mse;
      xhat = x;
      lambdastar = lambdas(i);
      epsilonstar = epsilons(j);
    end
  end
end
