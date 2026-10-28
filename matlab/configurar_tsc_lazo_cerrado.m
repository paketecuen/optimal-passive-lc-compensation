function cfg = configurar_tsc_lazo_cerrado(opciones)
%CONFIGURAR_TSC_LAZO_CERRADO Parametros y excitacion, sin simular maniobras.
% No llama al motor MATLAB ni genera estados/ordenes de control futuros.
% La carga es una inyeccion sintetica exogena; no tiene estados internos.
arguments
    opciones.CasoTension (1,1) string {mustBeMember(opciones.CasoTension,["original","moderada"])} = "moderada"
    opciones.Modo (1,1) string {mustBeMember(opciones.Modo,["ga","convencional","fijo"])} = "ga"
    opciones.BitsFijos (1,2) double = [1 0]
    opciones.PNominal_W (1,1) double {mustBePositive,mustBeFinite} = 1000
    opciones.Duracion_s (1,1) double {mustBePositive,mustBeFinite} = .8
    opciones.PasoMaximo_s (1,1) double {mustBePositive,mustBeFinite} = 2.5e-5
    opciones.VCInicial_V (1,2) double {mustBeFinite} = [0 0]
    opciones.MaxDiDt_A_s (1,1) double {mustBeNonnegative,mustBeFinite} = 2500
    opciones.Permanencia_s (1,1) double {mustBeNonnegative,mustBeFinite} = .06
    opciones.Histeresis_A2 (1,1) double {mustBeNonnegative,mustBeFinite} = .025
    opciones.MargenDiDt (1,1) double {mustBeNonnegative,mustBeFinite,mustBeLessThan(opciones.MargenDiDt,1)} = .01
    opciones.EscalaCarga (1,:) double {mustBePositive,mustBeFinite} = [.45 1.35 .55 1.15]
    opciones.InstantesCarga_s (1,:) double {mustBeNonnegative,mustBeFinite} = [0 .2 .4 .6]
    opciones.Transicion_s (1,1) double {mustBePositive,mustBeFinite} = .06
end
assert(all(ismember(opciones.BitsFijos,[0 1])),'TSCLazo:Bits','Se requieren bits binarios.');
assert(numel(opciones.EscalaCarga)==numel(opciones.InstantesCarga_s) && ...
    opciones.InstantesCarga_s(1)==0 && all(diff(opciones.InstantesCarga_s)>=opciones.Transicion_s), ...
    'TSCLazo:Carga','Las rampas deben estar ordenadas y no solaparse.');
c=caso_m3(63); scaled=dimensionar_referencia_230V(opciones.PNominal_W);
coef=[1 .03 .03]; if opciones.CasoTension=="original", coef=[1 1 1]; end
co=cos(c.tau); u0=sqrt(2)*(coef(1)*co+coef(2)*(2*co.^2-1)+coef(3)*(4*co.^3-3*co));
N=numel(co); primitive=[c.D,ones(N,1);ones(1,N),0]\[u0;0];
i0=c.G*u0+c.C_normalizada*c.D*u0+c.gamma_L*primitive(1:end-1)+ ...
    c.alpha_serie*((c.D*c.D+c.z_serie*eye(N))\(c.D*u0));
b=scaled.banco; weight=[.5 1]; T=2*pi/c.omega;
p=struct('L_H',b.Ls_H./weight,'C_F',b.Cs_F*weight, ...
    'R_ohm',c.omega*b.Ls_H./weight/100+1./(c.omega*b.Cs_F*weight*500), ...
    'Rb_ohm',1000*scaled.escala_impedancia./weight, ...
    'Cp_F',b.Cp_F,'Rp_ohm',1/(c.omega*b.Cp_F*500),'omega',c.omega,'T',T, ...
    'lambda',230/norm(coef),'k',scaled.escala_impedancia,'coefTension',coef, ...
    'limiteDiDt_A_s',opciones.MaxDiDt_A_s, ...
    'limiteDiDtPlan_A_s',(1-opciones.MargenDiDt)*opciones.MaxDiDt_A_s, ...
    'margenDiDt',opciones.MargenDiDt,'permanenciaMin_s',opciones.Permanencia_s, ...
    'histeresis_A2',opciones.Histeresis_A2,'pasoMuestreo_s',T/N, ...
    'ruidoSigma_A',0,'PNominal_W',opciones.PNominal_W,'opciones',opciones);
p.uReferencia=u0*p.lambda; p.iCargaReferencia=i0*p.lambda/p.k;
Dt=p.omega*c.D; icp=(eye(N)+p.Rp_ohm*p.Cp_F*Dt)\(p.Cp_F*Dt*p.uReferencia);
p.cpReferencia=[icp,p.uReferencia-p.Rp_ohm*icp];
p.x0=[0;opciones.VCInicial_V(1);0;opciones.VCInicial_V(2);0;0];
p.bits=[0 0;0 1;1 0;1 1];
cfg=struct('p',p,'D63',c.D,'Modo',opciones.Modo,'BitsFijos',opciones.BitsFijos, ...
    'TsSupervision',T/N,'Duracion_s',opciones.Duracion_s, ...
    'PasoMaximo_s',opciones.PasoMaximo_s,'VCInicial_V',opciones.VCInicial_V, ...
    'CasoTension',opciones.CasoTension,'opciones',opciones);
% Solo entrada exogena: ninguna decision ni maniobra precalculada.
t=(0:min(1e-5,opciones.PasoMaximo_s):opciones.Duracion_s)';
if t(end)<opciones.Duracion_s, t(end+1)=opciones.Duracion_s; end
scale=repmat(opciones.EscalaCarga(1),size(t));
for j=2:numel(opciones.EscalaCarga)
    phase=max(0,min(1,(t-opciones.InstantesCarga_s(j))/opciones.Transicion_s));
    scale=scale+(opciones.EscalaCarga(j)-opciones.EscalaCarga(j-1))*(1-cos(pi*phase))/2;
end
cfg.CargaExogena=[t,scale.*interpolar_periodico(p.iCargaReferencia,p.omega*t)];
cfg.vCp0=p.cpReferencia(1,2);
cfg.alcance="Solo simulacion; fuente impuesta, carga sintetica, valvula idealizada, L/C fijos. Control calculado en Simulink durante la ejecucion.";
end
