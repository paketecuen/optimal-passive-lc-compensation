function q = auditar_arranque_fuente_tsc(r)
%AUDITAR_ARRANQUE_FUENTE_TSC Energia incremental del sistema acoplado.
% Auditoria posterior, no nuevo controlador. Solo desde el primer ON hasta
% la siguiente maniobra/rampa. Se compara con la orbita periodica de esa
% topologia y esa carga constante, incluyendo los estados de la fuente.
a=r.trazas; p=r.parametros; cfg=r.config; events=r.eventos;
first=find(events.Tipo=="ON",1);
assert(~isempty(first),'TSCFuente:Arranque','No hay un encendido que auditar.');
t0=events.Tiempo_s(first); later=events.Tiempo_s(events.Tiempo_s>t0+1e-10 & ...
    ismember(events.Tipo,["ON","OFF"]));
stop=min([cfg.Duracion_s;later;cfg.opciones.InstantesCarga_s(2:end)']);
assert(t0<stop,'TSCFuente:Arranque','El primer ON no pertenece al primer tramo constante.');
% El solver puede registrar el estado abierto unos ulps antes del evento.
% Se toma el primer estado efectivamente conectado, sin desplazar el ON.
start=find(a.t>=t0-1e-12 & a.estadoAplicado(:,events.Rama(first))==1,1);
indices=find((1:numel(a.t))'>=start & a.t<stop-1e-10); t=a.t(indices);
assert(abs(t(1)-t0)<1e-9,'TSCFuente:Arranque','Falta el estado del primer encendido.');
bits=a.estadoAplicado(indices(1),:);
assert(all(a.estadoAplicado(indices,:)==bits,'all'),'TSCFuente:Arranque','Cambia la topologia dentro del tramo.');
ref=resolver_periodico_fuente_tsc(cfg,cfg.opciones.EscalaCarga(1),bits);
% Cada rama tiene psi_j=sqrt(L_j)*i_j P+ + sqrt(C_j)*v_j N-.
% El estado completo es la suma directa de esos planos y los estados de
% Cp/Lg. La norma cuadrada del desajuste es dos veces su energia.
values=[ref.iFuente,ref.vCp,ref.iCp,ref.iRamas,ref.vCs];
rp=interpolar_periodico(values,p.omega*t);
diG=a.iRed(indices)-rp(:,1); dvp=a.vCp(indices)-rp(:,2);
dicp=a.iCp(indices)-rp(:,3); di=a.iRamas(indices,:)-rp(:,4:5);
dv=a.vCs(indices,:)-rp(:,6:7);
energy=.5*cfg.Fuente.L_H*diG.^2+.5*p.Cp_F*dvp.^2 ...
    +.5*sum(p.L_H.*di.^2+p.C_F.*dv.^2,2);
loss=cfg.Fuente.R_ohm*diG.^2+p.Rp_ohm*dicp.^2 ...
    +sum(p.R_ohm.*di.^2+dv.^2./p.Rb_ohm,2);
residual=energy-energy(1)+cumtrapz(t,loss);
vmax=zeros(1,2); imax=vmax;
for j=1:2
    vmax(j)=maxPeriodic(ref.vCs(:,j)); imax(j)=maxPeriodic(ref.iRamas(:,j));
end
vcBound=vmax+sqrt(2*energy(1)./p.C_F);
iBound=imax+sqrt(2*energy(1)./p.L_H);
peakVC=max(abs(a.vCs(indices,:)),[],1); peakI=max(abs(a.iRamas(indices,:)),[],1);
q=struct('t0_s',t0,'fin_exclusivo_s',stop,'bits',bits,'E_desajuste_inicial_J',energy(1), ...
    'Cota_VC_V',vcBound,'Pico_VC_observado_V',peakVC, ...
    'Cota_I_A',iBound,'Pico_I_observado_A',peakI, ...
    'Residuo_balance_incremental_J',max(abs(residual)), ...
    'Crecimiento_energia_J',max(energy-energy(1)), ...
    't',t,'energiaDesajuste_J',energy,'perdidasDesajuste_W',loss, ...
    'alcance',"Cota condicional para topologia y excitaciones periodicas fijas. Los maximos de la referencia se estiman numericamente. No es una proteccion de la trayectoria completa.");
assert(q.Residuo_balance_incremental_J<5e-4*max(1,energy(1)), ...
    'TSCFuente:EnergiaIncremental','No cierra la energia del desajuste acoplado.');
assert(all(peakVC<=vcBound+1e-3) && all(peakI<=iBound+1e-4), ...
    'TSCFuente:CotaArranque','La trayectoria supera la envolvente estimada.');
end

function peak=maxPeriodic(values)
% Maximo continuo buscado numericamente sobre el interpolante temporal.
phase=linspace(0,2*pi,1001)'; v=abs(interpolar_periodico(values,phase));
peak=max(v); extrema=find(v(2:end-1)>=v(1:end-2) & v(2:end-1)>=v(3:end))+1;
for k=extrema'
    [~,negative]=fminbnd(@(tau)-abs(interpolar_periodico(values,tau)), ...
        phase(k-1),phase(k+1),optimset('TolX',1e-12));
    peak=max(peak,-negative);
end
end
