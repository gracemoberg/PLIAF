%{ 
Passive Leaky Integrate and Fire Model

T_m dV/dt = e_L - V+R_mI_e

T_m = 10 ms 
    - tau_m, a time constant (r_m*c_m) 
    - assumed to be positive 
    - a property of the cell
V(t): voltage
V_th = -50 mV
    - voltage threshold for a spike
e_L = V_reset = -65 mV
    - equilibrium w/respect to leak
I_e: input current, mA
    - typically not dependent on voltage
R_m = 10 Mohms
    - total membrane resistance 
    - measure of how difficult it is for I_e to change voltage
    - a property of the cell

***values of constants taken from Dayan Abbott 5.4 p12-13
%}

% set parameter values
par.Tm = 10; % membrane time constant
par.eL = -65; % equilibrium level
par.Rm = 10e6; % total membrane resistance
par.Ie = 10e-6; % constant input current
par.V_th = -50; % threshold voltage
par.V_spike = 0; % spike voltage
par.r_m = 1;
par.Tsra = 100;
par.Deltag = 0.5;
par.Ek = -90;

% initialize time vec
h = 0.25; % time step
tfinal = 100; % total length of time in ms
t = 0:h:tfinal;
iter = length(t);

% solns: first row voltage, second row is spike rate g_sra
X = zeros(2, iter); 

% initial condition values
X(1,1) = -55; % voltage initial condition
X(2,1) = 1; %g_sra initial condition


% Runge-Kutta 4
for k = 1:iter -1
    
    k1 = SRAPLIAF( X(:,k), par);
    k2 = SRAPLIAF( X(:,k) + k1.*h./2, par);
    k3 = SRAPLIAF( X(:,k) + k2.*h./2, par); 
    k4 = SRAPLIAF( X(:,k) + k3.*h, par);
    
    X(:,k+1) = X(:,k) + h.*( k1 + 2.*k2 + 2.*k3 + k4)./6;
    
    if X(1, k+1) > par.V_th
        X(1, k) = par.V_spike; % insert spike
        X(1, k+1) = par.eL; % reset to equilibrium level
        X(2, k+1) = X(2, k) + par.Deltag;
    end    
      
end

figure(1); hold on;
plot(t, X(1,:), 'linewidth', 2);
xlabel('Time (ms)'); ylabel('Voltage (ms)')
title('Leaky Integrate and Fire Model with Spike Rate Adaptation');
set(gca, 'fontsize', 18, 'linewidth', 2); box on;
xlim([0, tfinal]);


figure(2); hold on;
plot(t, X(2,:), 'linewidth', 2);
xlabel('Time (ms)'); ylabel('Spike Rate')
title('Spike Rate Adaptation');
set(gca, 'fontsize', 18, 'linewidth', 2); box on;
xlim([0, tfinal]);








