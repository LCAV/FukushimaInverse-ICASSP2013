function S = reconstructSourceL1Pos(A, b, lambda)
% Solve the L1 regularized problem with positivity constraint
%
% min |Ax-b|^2 - lambda |x|_1
% subject to x >= 0
%
% A : transport matrix
% b : measurement vector
% lambda : regularization parameters
%
% Return
% S : solution of the optimization problem


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


n = size(A, 2);

cvx_begin; 

  variable xx(n);

  minimize ((A*xx-b)'*(A*xx-b)+lambda*norm(xx,1)); %l1-norm minimization

  subject to
    xx >= 0 ;

cvx_end; 

% Solution
S = xx;

