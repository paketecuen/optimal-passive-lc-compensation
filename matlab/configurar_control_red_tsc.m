function cfg = configurar_control_red_tsc(cfg,opciones)
%CONFIGURAR_CONTROL_RED_TSC Predictor temporal acoplado y admision de ON.
% Requiere e medida aguas arriba, u/is/corrientes de rama y estados actuales.
% Limites de ensayo; no corresponden a una ficha de componente comercial.
arguments
    cfg (1,1) struct
    opciones.MaxVC_V (1,2) double {mustBePositive} = [1100 1100]
    opciones.MaxI_A (1,2) double {mustBePositive} = [8 8]
    opciones.ModeloCoste (1,1) string {mustBeMember(opciones.ModeloCoste,["red","bornes"])} = "red"
    opciones.Margen (1,1) double {mustBeNonnegative,mustBeLessThan(opciones.Margen,1)} = .01
    opciones.PasosPlan (1,1) double {mustBeInteger,mustBePositive} = 126
    opciones.MuestrasCota (1,1) double {mustBeInteger,mustBePositive} = 1008
end
assert(isfield(cfg,'Fuente'),'TSCRed:Fuente','Configure primero la fuente R-L.');
assert(cfg.Modo=="ga",'TSCRed:Modo','Esta extension usa el criterio RMS temporal.');
cfg.NombreModelo="tsc_control_red";
c=opciones; p=cfg.p; N=63; I=eye(N); D=p.omega*cfg.D63;
c.Z=cfg.Fuente.R_ohm*I+cfg.Fuente.L_H*D;
c.Ycp=(I+p.Rp_ohm*p.Cp_F*D)\(p.Cp_F*D);
c.Yi=cell(1,2); c.Yv=cell(1,2);
for j=1:2
    H=[p.L_H(j)*D+p.R_ohm(j)*I,I;-I,p.C_F(j)*D+I/p.Rb_ohm(j)]\[I;zeros(N)];
    c.Yi{j}=H(1:N,:); c.Yv{j}=H(N+1:end,:);
end
c.pesos=[cfg.Fuente.L_H,p.Cp_F,p.L_H(1),p.C_F(1),p.L_H(2),p.C_F(2)];
c.escala=sqrt(c.pesos); if c.escala(1)==0, c.escala(1)=1; end
c.retardos=(0:opciones.PasosPlan)'*p.T/opciones.PasosPlan;
c.Kplan=interpolar_periodico(I,p.omega*c.retardos);
c.Kcota=interpolar_periodico(I,(0:opciones.MuestrasCota-1)'*2*pi/opciones.MuestrasCota);
c.topologias=cell(4,1);
for k=1:4
    bits=p.bits(k,:); Y=c.Ycp+bits(1)*c.Yi{1}+bits(2)*c.Yi{2};
    A=modelo_estado(cfg,bits); An=(c.escala'.*A)./c.escala;
    F=zeros(6,6,numel(c.retardos));
    for j=1:numel(c.retardos), F(:,:,j)=expm(An*c.retardos(j)); end
    c.topologias{k}=struct('Ue',(I+c.Z*Y)\I,'Ul',-(I+c.Z*Y)\c.Z, ...
        'Y',Y,'A',A,'F',F);
end
cfg.ControlRed=c;
cfg.alcance="Simulacion con predictor temporal de red y limites condicionales de arranque. e se mide aguas arriba; carga reconstruida localmente. Se repite el ultimo periodo, sin anticipar rampas. No garantiza limites frente a excitaciones futuras arbitrarias.";
end

function A=modelo_estado(cfg,bits)
% x=[is,vCp,i1,vC1,i2,vC2]. Para Lg=0, x1 es una coordenada nula.
p=cfg.p; f=cfg.Fuente; I=eye(6); A=zeros(6);
if f.L_H>0
    U=[p.Rp_ohm,1,-p.Rp_ohm,0,-p.Rp_ohm,0];
    A(1,:)=(-f.R_ohm*I(1,:)-U)/f.L_H;
else
    U=[0,f.R_ohm,-f.R_ohm*p.Rp_ohm,0,-f.R_ohm*p.Rp_ohm,0]/(p.Rp_ohm+f.R_ohm);
end
A(2,:)=(U-I(2,:))/(p.Rp_ohm*p.Cp_F);
for j=1:2
    ii=2*j+1; vv=ii+1;
    A(ii,:)=bits(j)*(U-p.R_ohm(j)*I(ii,:)-I(vv,:))/p.L_H(j);
    A(vv,:)=(bits(j)*I(ii,:)-I(vv,:)/p.Rb_ohm(j))/p.C_F(j);
end
end
