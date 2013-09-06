function [xhat, err, lambdastar] = tikhonov(A, bn, s, lambdas)
% Solves the inverse problem using Tikhonov regularization
%     min |A*x - bn|^2 + lambda*|x|^2
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


n = length(s);

for i=1:length(lambdas)
  %Analytic solution
  x=(A'*A+lambdas(i)*eye(n))\(A'*bn);
  
  % Check for minimum over lambdas
  mse=mean((x-s).^2);
  if (i == 1 | mse < err)
    err = mse;
    xhat = x;
    lambdastar = lambdas(i);
  end
end
