function [sys,x0,str,ts]=RST_switch(t,x,u,flag,R1,S1,R2,S2,R3,S3,Ts,nr,ns,nt)
% This S-Function computes the control signal of multi RST controllers.
% There are three fixed RST controllers corresponding to the switching
% signal {1,2,3} and one adaptive RST controller (switch=4).
% The input vector u has three input signals for the switching signal,
% reference signal and the plant output and three vector inputs for R, S
% and T polynomials (that are computed from an adaptive model outside the block).
% The output of the block is the control signal. Ts is the sampling period.


% It is assumed that all controller have the same order, otherwise the code
% should be modified.

n=nr+ns+nt;

%u(1): Switching signal
%u(2): reference signal
%u(3): output signal

persistent y_old u_old

%Usatplus=2;Usatminus=-2;

switch flag,
    % Initialization
    case 0,
        sizes = simsizes;
        sizes.NumContStates  = 0;
        sizes.NumDiscStates  = n;
        sizes.NumOutputs     = 1;
        sizes.NumInputs      = n+6;
        sizes.DirFeedthrough = 1;
        sizes.NumSampleTimes = 1;
        
        sys = simsizes(sizes);
        
        x0  = zeros(sizes.NumDiscStates,1);
        str = [];
        ts  = [Ts 0];
        y_old = zeros(1,max([length(R1),length(R2),length(R3)]));
        u_old = zeros(1,max([length(S1),length(S2),length(S3)])-1);
        % state update
    case 2
        switch u(1)
            case 1
                R=R1;
                S=S1;
                % RST controller for G1
                T=sum(R1);
            case 2
                R=R2;
                S=S2;
                % RST controller for G2
                T=sum(R2);
            case 3
                R=R3;
                S=S3;
                % RST controller for G3
                T=sum(R3);
            case 4
                R=u(4:nr+4)';
                S=u(nr+5:nr+ns+5)';
                T=u(nr+ns+6:nr+ns+nt+6)';
        end
        % Compute the control signal and save it (do not forget the saturation)
        ind_r = 1:length(R);
        ind_s = 1:length(S)-1;
        
        r_old = u(2);
        y_old = [u(3),y_old(1:end-1)];
        temp = addPoli(T.*r_old,R.*(-y_old(ind_r)));
        u_k = sum(addPoli(temp,S(2:end).*(-u_old(ind_s))))/S(1);
        u_old = [u_k,u_old(1:end-1)];
        % Update the state vector (including past inputs, past outputs and past reference signals)
        if nt>0
            sys=[u_k;x(1:ns-1);u(3);x(ns+1:nr+ns-1);u(2);x(nr+ns+1:n-1)];
        else
            sys=[u_k;x(1:ns-1);u(3);x(ns+1:nr+ns-1)];
        end
        % output update
    case 3
        u_k = x(1);
        % Compute again u_k and send it out
        if(sum(isnan(u_k)))
            a =1;
        end
        sys=u_k;
        
    case 9
        sys=[];
end


