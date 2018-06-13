%--------------------------------------------------------------------------
%% Function Name: poleplace
%
% Assumptions: None
%
% Inputs:
%   A: vector with coefficients of polynomial A
%   B: vector with coefficients of polynomial B with delay
%   Hr: vector with coefficients of polynomial Hs
%   Hs: vector with coefficients of polynomial Hr
%   P: vector with coefficients of polynomial of desired poles P
% Outputs:
%   R: row vector of coefficients R
%   S: row vector of coefficients S

% $Date: April 25, 2018
%--------------------------------------------------------------------------
function [R,S]=poleplace(B,A,Hr,Hs,P)
% Compute degree of S', R' and P
    A = row2col(A);B = row2col(B);Hr = row2col(Hr);Hs = row2col(Hs);
    % B already contains the delay
    A1 = conv(A,Hs); 
    B1 = conv(B,Hr);
    nR1 = degree(A1) - 1;
    nS1 = degree(B1) - 1;
% Compute silvester Matrix
     M_A = zeros(nS1+1+nR1+1,nS1+1);
     M_B = zeros(nS1+1+nR1+1,nR1+1);

%     M_A = zeros(degree(A1) + degree(B1),degree(B1));
%     M_B = zeros(degree(A1) + degree(B1),degree(A1));
    for j = 1: nS1+1
        M_A(j:degree(A1)+j,j) = A1;
    end
    for j = 1: nR1+1
        M_B(end-j-degree(B1)+1:end-j+1,end-j+1) = B1;
    end
    M = [M_A,M_B];
% Compute R and S coefficients
    p = zeros(degree(A1) + degree(B1),1);
    P = row2col(P);
    p(1:length(P),1) = P;
    x = M\p;
    S = conv(Hs,x(1:nS1+1))';
    R = conv(Hr,x(nS1+2:end))';
end
function [deg] = degree(P)
    deg = length(P)-1;
end
function [A] = row2col(A)
if(isrow(A)), A = A'; end;
end