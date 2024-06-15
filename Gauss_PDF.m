% GAUSS_PDF  Multivariate Gaussian PDF
%
%  Syntax:
%   [P,E] = GAUSS_PDF(X,M,S)
%
%  In:
%   X - Dx1 value or N values as DxN matrix
%   M - Dx1 mean of distibution or N values as DxN matrix.
%   S - DxD covariance matrix
%
%  Out:
%   P - Probability of X. 
%   E - Negative logarithm of P
%   
%  Description:
%   Calculate values of PDF (Probability Density
%   Function) of multivariate Gaussian distribution
%


function [P,E] = Gauss_PDF(X,M,S)

  if size(M,2) == 1
    DX = X-repmat(M,1,size(X,2));
    E = 0.5*sum(DX.*(S\DX),1);  
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  elseif size(X,2) == 1
    DX = repmat(X,1,size(M,2))-M;  
    E = 0.5*sum(DX.*(S\DX),1);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  else
    DX = X-M;  
    E = 0.5*DX'*(S\DX);
    d = size(M,1);
    E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
    P = exp(-E);
  end

end