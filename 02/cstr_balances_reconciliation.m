function [c,ceq] = cstr_balances_reconciliation(y,V,CAin,Tin,Tc,...
                                                rho,Cp,k0,E,R,dH,UA)

Fin  = y(1);
Fout = y(2);
CA   = y(3);
T    = y(4);

k = k0*exp(-E/(R*T));
rA = k*CA;

%% 1. Balanço global de massa
eq_mass = Fin - Fout;

%% 2. Balanço de componente A
% eq_comp = Fin*CAin - Fout*CA - V*rA;
eq_comp = Fin*(CAin - CA) - V*rA;

%% 3. Balanço de energia
eq_energy = Fin*rho*Cp*(Tin-T)+ (-dH)*V*rA - UA*(T - Tc);

ceq = [eq_mass;
       eq_comp;
       eq_energy];

c = [];

end