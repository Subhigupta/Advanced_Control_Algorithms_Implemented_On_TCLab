function dDt = CSTRModel3(t, x, u)

Ca = x(1);
T = x(2);

Fl = u(1);
Fc = u(2);

% Parameters
kr= 6;%l/h
Vl=24;%l
rhol=800;%kg/m3
rhoc=1000;%kg/m3
Cpl=3;%kj/kg-K
Cpc=4.19;
Tl0=283;
Tc0=273;
Tc=373;
CA0=4;

%% Equations
dDdt(1,1) = (Fl*rhol*Cpl*(Tl0-Tl) + Fc*rhoc*Cpc*(Tc0 - Tc)+ Vl*kr*(CA0-Cb)*H) / (Vl*rho_l*Cp_l);
dDdt(2,1) = (Vl*kr*(CA0-Cb) - Fl*Cb) / Vl;