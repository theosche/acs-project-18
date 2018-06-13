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
        
        x0  = [zeros(m,1);1;DT];
        str = [];
        ts  = [-1 0]; 
        
     % state update 
    case 2
        % Based on u(i) compute the monitoring signal Ji in a recursive
        % way for all inputs (prediction errors)
        
        for i=1:m
            Ji_longterm = u(i).^2 + exp(-lambda)*x(1:m);
            Ji = u(i).^2 + beta * Ji_longterm;
        end
       
        % Write an algorithm for the Dwell-Time
        % DT shows the number of sampling period of waiting between two
        % switchings and is saved in x(m+2).
        % x(m+1) is the choice of the best predictor.
        if (x(m+2) < DT)
            x(m+2) = x(m+2) + 1; % Increase current time
        end
        
        [~, best] = min(Ji); % Instantaneous best controller
        if (x(m+2) == DT && best ~= x(m+1))
            x(m+2) = 0;
            x(m+1) = best;
        end
        
        sys=[Ji_longterm;x(m+1);x(m+2)];
        
     % output update
    case 3
        sys=x(m+1); 
    case 9
        sys=[];
 end
    
        
        