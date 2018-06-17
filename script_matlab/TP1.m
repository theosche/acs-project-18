%% TP1 Hinf Modifie
close all;
z = tf('z',Ts); % Set z as variable
s = tf('s'); % Set s as variable

% Parameters
n = 4; % Order of the basis -> to change
a = 0; % Understood parameter -> to change
wc = 3; % Define bandwidth of our system -> to change
wh = 50;
Mm = 0.7; % Modulus Margin

Ld = wc/s; % Desired L
phi = conphi('Laguerre',[Ts a n],'z',z/(z-1)); % Laguerre basis
per = conper('LS',Mm,Ld); % Set perfomance -> minimize 2-norm(L-Ld) with Mm
K3 = condes(G1, phi, per);
% NyquistConstr(K3, G1, per); axis square % Evaluate modulus margin constraint
% T3 = feedback(G1*K3,1); % Evaluate settling time
% temp = stepinfo(T3);
% tau= temp.SettlingTime;
% disp(tau);
% figure(); step(T3); grid on
% figure(); bode(G1*K3);grid on