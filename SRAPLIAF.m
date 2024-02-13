function Xdot = SRAPLIAF(X, par)

V = X(1);
g = X(2);

fV = (1/par.Tm) * (par.eL - V - par.r_m*g*(V - par.Ek) + par.Rm*par.Ie);
fg = (1/par.Tsra) * (-g);

Xdot = [fV; fg];