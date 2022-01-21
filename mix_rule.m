function f = mix_rule(epsh,epsr, epsf, nu)
%MIX_RULE calculate volume fraction based on mixing rule
%   with calculated dielectric constant 
%epsh = dielectric of half-space (background)
%epsr = dielectric of rock (inclusions)
%epsf = measured effective dielectric constant
%nu = mixing model identifier
%   0 = Maxwell-Garnet 
%   2 = Polder van-Santen
%   3 = Coherent potential approx. 
%%%%%epsi is ice, could possibly change a lot with temperature.

%%%Reference: Shivola, 1999: EM Mixing Formulas & Applications, section
%%%9.3.1

A = (epsf-epsh)./(epsf+2*epsh+nu.*(epsf-epsh));
B = (epsr-epsh)./(epsr+2*epsh+nu.*(epsf-epsh));
f = A./B;

end

