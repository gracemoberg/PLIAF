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

% set parameter values
Tm = 10; % membrane time constant
eL = -65; % equilibrium level
Rm = 10e6; % total membrane resistance
Ie = 1.6e-6; % constant input current
V_th = -50; % threshold voltage
V_spike = 0; % spike voltage

% set to 0 for constant input current, 1 for piecewise
pwCurrent = 0; % flag

% if defining a piecewise input current
Ie1 = 0; % first quarter input current
Ie2 = 1.5e-6;% second quarter input current
Ie3 = 5e-6; % third quarter input current
Ie4 = 1.7e-6; % fourth quarter input current

% initial voltage
V0 = -55; %initial condition

% initialize time vec
h = 0.25; % time step
tfinal = 100; % total length of time
t = 0:h:tfinal;

% construct basic function f
f = @(t,V)(1/Tm*(eL - V + Rm*Ie));

% initialize numerical approximation
soln = zeros(1,length(t));
soln(1) = V0;

% simulate using Runge-Kutta 4

if pwCurrent == 1 % if piecewise flag = 1
    for n = 1:length(t)-1
        
        if n < (length(t)-1)/4 % first quarter
            Ie = Ie1; 
            f = @(t,V)(1/Tm*(eL - V + Rm*Ie));
        elseif (n >= (length(t)-1)/4) && (n < (length(t)-1)/2) % second quarter
            Ie = Ie2; 
            f = @(t,V)(1/Tm*(eL - V + Rm*Ie));
        elseif (n >= (length(t)-1)/2) && (n < 3*(length(t)-1)/4) % third quarter
            Ie = Ie3; 
            f = @(t,V)(1/Tm*(eL - V + Rm*Ie));
        else % fourth quarter
            Ie = Ie4; 
            f = @(t,V)(1/Tm*(eL - V + Rm*Ie));
        end

        k1 = f(t(n), soln(n));
        k2 = f(t(n)+h/2, soln(n)+ h*k1/2);
        k3 = f(t(n)+h/2, soln(n)+ h*k2/2);
        k4 = f(t(n)+h, soln(n)+ h*k3);

        % Solution update
        soln(n + 1) = soln(n) + h*( k1 + 2*k2 + 2*k3 + k4)/6;

        % insert spike if exceeds threshold
        if soln(n+1) > V_th
            soln(n) = V_spike; % insert spike
            soln(n+1) = eL; % reset to equilibrium level
        end    
    end
    
else 
    for n = 1:length(t)-1

        k1 = f(t(n), soln(n));
        k2 = f(t(n)+h/2, soln(n)+ h*k1/2);
        k3 = f(t(n)+h/2, soln(n)+ h*k2/2);
        k4 = f(t(n)+h, soln(n)+ h*k3);

        % Solution update
        soln(n + 1) = soln(n) + h*( k1 + 2*k2 + 2*k3 + k4)/6;

        % insert spike if exceeds threshold
        if soln(n+1) > V_th
            soln(n) = V_spike; % insert spike
            soln(n+1) = eL; % reset to equilibrium level
        end
        
    end
end
    
plot(t, soln)
title('Leaky Integrate and Fire Model')
xlabel('t (ms)')
ylabel('V (mV)')









