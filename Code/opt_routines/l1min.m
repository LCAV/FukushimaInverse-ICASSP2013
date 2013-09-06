function [xhat, err, lambdastar] = l1min(A, bn, s, lambdas, scaling)
% Solves the inverse problem using Tikhonov regularization
%     min |A*x - bn|_2^2 + lambda*|x|_1
% for all given values of lambda.
%
% Param:
%   A       : forward matrix
%   bn      : noisy measurement vector
%   s       : original source vector (bn = A*s+noise)
%   lambdas : values of the optimization parameter
%   scaling : a pre-conditioning parameter to ensure the stability of the problem
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

A = scaling*A;

n=size(s,1); % length of the source vector

for i=1:size(lambdas,2)
  % Solve optimization problem with CVX
  cvx_begin; 
  cvx_precision high
    variable x(n);
    minimize((A*x-bn)'*(A*x-bn)+lambdas(i)*norm(x,1)); %l1-norm minimization
  cvx_end;  

  % Check for minimum over lambdas
  mse=mean((scaling*x-s).^2);
  if (i == 1 || mse < err)
    err = mse;
    xhat = x;
    lambdastar = lambdas(i);
  end
end
