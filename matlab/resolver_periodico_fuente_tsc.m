function r = resolver_periodico_fuente_tsc(cfg,escala,bits)
%RESOLVER_PERIODICO_FUENTE_TSC Auditoria estacionaria mediante operadores temporales.
% En el espacio de muestras, D representa el derivador antiautoadjunto.
% u=(I+Zg*Yb)^-1(e-Zg*iCarga), is=iCarga+Yb*u. No usa FFT.
% Solo respuesta periodica forzada para una topologia y carga constantes.
arguments
    cfg (1,1) struct
    escala (1,1) double {mustBePositive,mustBeFinite}
    bits (1,2) double
end
assert(all(ismember(bits,[0 1])),'TSCFuente:Bits','Se requieren bits binarios.');
p=cfg.p; N=size(cfg.D63,1); I=eye(N); D=p.omega*cfg.D63;
Z=cfg.Fuente.R_ohm*I+cfg.Fuente.L_H*D;
Ycp=(I+p.Rp_ohm*p.Cp_F*D)\(p.Cp_F*D); Y=Ycp;
Yi=cell(1,2); Yv=cell(1,2);
for j=1:2
    H=[p.L_H(j)*D+p.R_ohm(j)*I,I;-I,p.C_F(j)*D+I/p.Rb_ohm(j)]\[I;zeros(N)];
    Yi{j}=H(1:N,:); Yv{j}=H(N+1:end,:);
    Y=Y+bits(j)*Yi{j};
end
e=p.uReferencia; il=escala*p.iCargaReferencia;
u=(I+Z*Y)\(e-Z*il); icp=Ycp*u; is=il+Y*u;
ir=zeros(N,2); vc=ir;
for j=1:2, ir(:,j)=bits(j)*Yi{j}*u; vc(:,j)=bits(j)*Yv{j}*u; end
loss=p.Rp_ohm*icp.^2+sum(p.R_ohm.*ir.^2+vc.^2./p.Rb_ohm,2);
r=struct('u',u,'e',e,'iFuente',is,'iCarga',il,'iCp',icp, ...
    'vCp',u-p.Rp_ohm*icp,'iRamas',ir,'vCs',vc, ...
    'I_fuente_RMS_A',sqrt(mean(is.^2)),'U_RMS_V',sqrt(mean(u.^2)), ...
    'P_carga_W',mean(u.*il),'P_fuente_W',mean(e.*is), ...
    'P_perdidas_W',mean(loss)+cfg.Fuente.R_ohm*mean(is.^2), ...
    'Residuo_KVL_V',max(abs(e-u-Z*is)), ...
    'Residuo_KCL_A',max(abs(is-il-icp-sum(ir,2))), ...
    'Balance_W',mean(e.*is-u.*il-loss)-cfg.Fuente.R_ohm*mean(is.^2));
end
