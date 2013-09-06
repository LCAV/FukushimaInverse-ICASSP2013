function [P,V] = matrixCleaning(M, V, CN, lb, ub)
% [P,th] = MCleaning(M, CN, lb, ub)
% 
% M: matrix to clean
% V: associated measurement vector
% CN: desired condition number
% lb: lower bound for search  (optional)
% ub: upper bound for search  (optional)
%
% Return:
% P: cleaned M
% th: the threshold corresponding to CN


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


% these bounds will find the right M
if (nargin < 4)
  lb = 0.01;
  ub = 0.1;
elseif (nargin < 5)
  ub = 0.1;
end

% number of points to test
Npts = 200;

% take the thresholds to test equally spaced
ths = linspace(lb,ub,Npts);   %threasolds for the norm

% measure norms of all columns
m = sqrt(sum(M.*M, 2));
maxNorm = max(m);

% for every threshold tested, remove rows with (norm/maxNorm < threshold)
% stop when condition number is <= CN
P = M;
for thsIndx=1:size(ths,2)
  ind_2 = find(m > ths(thsIndx)*maxNorm);
  V = V(ind_2);
  P = P(ind_2,:);
  m = m(ind_2);
  cn = cond(P);
  %fprintf('thsIndx=%d th=%f cn=%f Nrow=%d\n', thsIndx, ths(thsIndx), cn, size(P, 1));
  if (cn <= CN)
    th = ths(thsIndx);
    break;
  end
end

