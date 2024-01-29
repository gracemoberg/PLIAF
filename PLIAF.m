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

***values of constants taken from Dayan Abbott 5.4 p12
%}

%{
V(t + dt) = V_inf + (V(t) - V_inf)*exp(-dt/tau_v)
%}


%  set parameter values
Tm = 10; 
eL = -65;
Rm = 10e6;
Ie = 2.2e-6;
V_th = -50;
V_spike = -20;


% initialize time vec
h = 0.25; 
tfinal = 100; 
t = 0:h:tfinal;
% squeeze(t);

% construct function f
f = @(t,V)(1/Tm*(eL - V + Rm*Ie));
soln0 = -55; %initial condition

% numerical approximation
soln = zeros(1,length(t));
size(soln)
soln(1) = soln0;

% size(soln)

% simulate using Runge-Kutta 4
for n = 1:length(t)-1
    k1 = f(t(n), soln(n));
    k2 = f(t(n)+h/2, soln(n)+ h*k1/2);
    k3 = f(t(n)+h/2, soln(n)+ h*k2/2);
    k4 = f(t(n)+h, soln(n)+ h*k3);
    
%   Solution update
    soln(n + 1) = soln(n) + h*( k1 + 2*k2 + 2*k3 + k4)/6;
    
    % insert spike if exceeds threshold
    if soln(n+1) > V_th
        soln(n) = V_spike;
        soln(n+1) = eL;
    end
        

end

size(t)
size(soln)

plot(t, soln)








