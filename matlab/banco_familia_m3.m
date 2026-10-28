function b = banco_familia_m3(f,z)
%BANCO_FAMILIA_M3 Un miembro fisico de la familia optima; ceros omiten ramas.
arguments
    f (1,1) struct
    z (1,1) double {mustBePositive,mustBeFinite}
end
assert(min(abs(z-f.modelo.polos))>1e-9*max(z,1), ...
    'LCfamilia:Resonancia','Resonancia excitada: no hay miembro finito.');
coef=[polyval(fliplr(f.c_numerador),z);polyval(fliplr(f.gamma_numerador),z)/z; ...
    polyval(fliplr(f.alpha_numerador),z)/z];
tol=1e-10*max(norm(coef),1);
assert(all(coef(1:2)>=-tol) && coef(3)>tol,'LCfamilia:Signos','Banco no pasivo en esta carta.');
coef(abs(coef)<tol)=0; % Solo redondeo de una frontera; no proyeccion por recorte.
Lp=Inf; if coef(2)>0, Lp=1/(f.omega*coef(2)); end
b=struct('coef',coef,'z',z,'Cp_F',coef(1)/f.omega,'Lp_H',Lp, ...
    'Ls_H',1/(f.omega*coef(3)),'Cs_F',coef(3)/(f.omega*z), ...
    'f_serie_Hz',f.omega*sqrt(z)/(2*pi),'componentes',sum(coef(1:2)>0)+2);
end
