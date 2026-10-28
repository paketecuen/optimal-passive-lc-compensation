function out = escalar_banco_m3(b,scale)
%ESCALAR_BANCO_M3 Otra red pasiva fija: misma sintonia y admitancia escalada.
arguments
    b (1,1) struct
    scale (1,1) double {mustBePositive,mustBeFinite}
end
out=b; out.coef=b.coef*scale;
if isfield(b,'Cp_F')
    out.Cp_F=b.Cp_F*scale; out.Lp_H=b.Lp_H/scale;
    out.Ls_H=b.Ls_H/scale; out.Cs_F=b.Cs_F*scale;
end
end
