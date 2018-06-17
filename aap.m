function [sys,x0,str,ts]=aap(t,x,u,flag,nA,nB,d)
% This function is a parametric adaptation algorithm. It takes the input
% and output of a system and computes the model parameters. The adaptation
% gain is initialized to F0 and the parameters are initialized to the
% parameters of the last model chosen by the switching algorithm. The
% function has three inputs: u(1) is the plant input, u(2) the plant output
% and u(3) is the switching signal belongs to {1,2,3,4}. The output of the
% block is the predicted output yhat and the vector of estimated
% parameters theta. The parameter vector theta is initialized by the best
% model among the fixed one after each transition of the switching signal u(3).
persistent sigma u_delay
n=nA+nB;
Nstate=n*(n+2);
deadzone=0.01; % this value can be chosen to freeze adaptation if the adaptation error is too small.
F0=eye(n);
alpha_gain = 1;
GI = 0.2;

% The parameters of G1
B1=[0         0   -0.0062    0.0118   -0.0052    0.0068]';
A1=[1.0000   -4.2442    8.4788  -10.3763    8.2172   -3.9835    0.9080]';

% The parameters of G2
B2=[0         0    0.0179   -0.0410    0.0475   -0.0074]';
A2=[1.0000   -3.9603    7.5101   -8.8773    6.9061   -3.3876    0.8090]';

% The parameters of G3
B3=[0         0    0.0080   -0.0118    0.0204    0.0079]';
A3=[1.0000   -3.7349    6.9589   -8.3407    6.6806   -3.4386    0.8750]';


switch flag,
    case 0,
        
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = Nstate;
        sizes.NumOutputs     = n+1;
        sizes.NumInputs      = 3;
        sizes.DirFeedthrough = 0;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        sigma = 0;
        u_delay = zeros(d+1,1);
        x0  = zeros(Nstate,1);
        str = [];
        ts  = [-1 0];
        
        theta_k = zeros(size([A1(2:end);B1(d+2:end)]));
        phi_p = zeros(1,n);
        F0 = GI*diag(ones(n,1));
        sys=[theta_k;phi_p;reshape(F0,n*n,1)];
        
    case 2
        if(~sigma)
            sigma = mod(u(3)-1,3)+1;
        end
        if(u(3)~= sigma) % Check if the sigma change
            switch(sigma)
                case 1
                    theta_k = [A1(2:end);B1(d+2:end)];
                case 2
                    theta_k = [A2(2:end);B2(d+2:end)];
                case 3
                    theta_k = [A3(2:end);B3(d+2:end)];
                case 4
                    theta_k = x(1:n);
            end
            x(2*n+1:end) = zeros(size(x(2*n+1:end)));
        else
            theta_k = x(1:n);
        end
        sigma = u(3);
        
        phi_k=[x(n+1:n+nA);x(n+nA+1:2*n)];
        F_k=reshape(x(2*n+1:end),n,n);
        
        % Compute F(t+1) using matrix inversion lemma
        
        %F_p = F_k - F_k*(phi_k*phi_k')*F_k/(1+phi_k'*F_k*phi_k);
        % With a constant trace adaptation gain
        F_p = F_k - F_k*(phi_k*phi_k')*F_k/(alpha_gain+phi_k'*F_k*phi_k);
        F_p = n*GI/trace(F_p)*F_p;
        
        % Compute the a priori prediction error
        epsilon = u(2) - theta_k'*phi_k;
        
        % Dead zone
        if epsilon < deadzone
            epsilon = 0;
        end
        
        % Compute the new estimate for theta
        theta_p = theta_k + F_p*phi_k*epsilon;
        
        % Update the observation vector
        u_delay = [u_delay(2:end); u(1)];
        phi_p=[-u(2);x(n+1:n+nA-1);u_delay(1);x(n+nA+1:2*n-1)];
        
        % Initialize the parameters with the best model in the transition
        % of the switching signal (its past value is saved in x(end))
        
 
        sys=[theta_p;phi_p;reshape(F_p,n*n,1)];
        
    case 3
        % Compute yhat and theta_k
        theta_k = x(1:n);
        phi_k = [x(n+1:n+nA);x(n+nA+1:2*n)];
        yhat = theta_k'*phi_k;
        sys=[yhat;theta_k];
        
    case 9
        sys=[];
end


