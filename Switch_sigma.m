function [sys,x0,str,ts]=Switch_sigma(t,x,u,flag,beta,lambda,DT)
% This S-Function computes the switching signal for switching adaptive control. 
% The input vector u has m inputs that are the estimation errors from a 
% multi-estimator. The parameters are the weighting factor beta 
% (alpha is taken equal to 1), forgetting factor lambda and Dwell-Time DT.


m=4;
x_p=zeros(m,1);

switch flag,
    % Initialization
    case 0,
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = m+2;
        sizes.NumOutputs     = 1;
        sizes.NumInputs      = m;
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;

        sys = simsizes(sizes);
        
        x0  = [zeros(m,1);0;1];
        str = [];
        ts  = [-1 0]; 
        
     % state update 
    case 2
        % Based on u(i) compute the monitoring signal Ji in a recursive
        % way for all inputs (prediction errors)
        
        Ji_longterm = u.^2 + exp(-lambda)*x(1:m);
        Ji = u.^2 + beta * Ji_longterm;
        
        % Write an algorithm for the Dwell-Time
        % DT shows the number of sampling period of waiting between two
        % switchings and is saved in x(m+2).
        % x(m+1) is the choice of the best predictor.
        
        x(m+1) = x(m+1) + 1; % Increase the counter 
        
        [~, best] = min(Ji(1:3)); % Instantaneous best controller
        if (x(m+1) == DT) % Dwell-Time
            x(m+1) = 0;
            x(m+2) = best;
        end
        
        sys=[Ji_longterm;x(m+1);x(m+2)];
        
     % output update
    case 3
        sys=x(m+2); 
    case 9
        sys=[];
 end
    
        
        