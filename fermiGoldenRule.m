function matelemsq = fermiGoldenRule(Sys,mag_field,transition)
% fermiGoldenRule - output | <b|S x|a> |^2
%   matelemsq = fermiGoldenRule(Sys,Exp,transition)
% 
% Sys describes spin system (max 1 elec and 1 nuc)
% mag_field is a column vector specifying mag field in crystal frame
% transition is a row vector specifying levels transitioning between
%
% This function is intended as an alternative to the Intensity output from
% easyspin's resfields and resfreqs functions
%
% Function currently does not include any of the other factors in the Fermi
% golden rule. Output is unitless.
%

H = sham(Sys,mag_field);
[V,E] = eig(H);
if any(strcmp('Nucs',fieldnames(Sys)))
    Sx = sop(Sys,'xe'); % Sx on electron
else
    Sx = sop(Sys,'x');
end
matelem = V(:,transition(1))'*Sx*V(:,transition(2));

matelemsq = norm(matelem)^2;