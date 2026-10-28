function informe = verificar_tsc_lazo_cerrado(r, opciones)
%VERIFICAR_TSC_LAZO_CERRADO Auditoria independiente de una corrida Simulink.
% Reconstruye costes, leyes electricas y energia desde las trazas guardadas.
% No llama al controlador ni vuelve a simular. Una entrada externa de carga
% es admisible; un cronograma externo de estados aplicados no lo es.
arguments
    r (1,1) struct
    opciones.Exportar (1,1) logical = false
    opciones.CarpetaSalida (1,1) string = ""
    opciones.InspeccionarModelo (1,1) logical = true
    opciones.AuditarSupervisor (1,1) logical = true
end
assert(all(isfield(r,{'trazas','parametros','config','eventos', ...
    'decisiones','ventanas_decision','modelo'})), ...
    'TSCLazo:Contrato','Registro incompleto para la auditoria.');
a=r.trazas; p=r.parametros; t=a.t(:); n=numel(t);
assert(n>2 && all(diff(t)>0),'TSCLazo:Tiempo','Se requiere tiempo estrictamente creciente.');
assert(all(isfinite([t;a.u(:);a.iCarga(:);a.iFuente(:);a.iCp(:); ...
    a.vCp(:);a.iRamas(:);a.vCs(:);a.integralEntrada_J(:);a.integralPerdidas_J(:)])), ...
    'TSCLazo:Finito','Las trazas fisicas contienen valores no finitos.');
assert(size(a.iRamas,2)==2 && size(a.vCs,2)==2, ...
    'TSCLazo:Ramas','La auditoria requiere dos LC conmutables.');
assert(all(ismember(a.estadoAplicado(:),[0 1])) && ...
    all(ismember(a.estadoSolicitado(:),[0 1])) && all(ismember(a.pendiente(:),[-1 0 1])), ...
    'TSCLazo:Bits','Estados aplicados, solicitados o pendientes invalidos.');
assert(~any(a.pendiente(:)==1 & logical(a.estadoAplicado(:))) && ...
    ~any(a.pendiente(:)==-1 & ~logical(a.estadoAplicado(:))), ...
    'TSCLazo:Pendiente','El estado pendiente contradice la conduccion real.');
L=reshape(p.L_H,1,[]); C=reshape(p.C_F,1,[]);
R=reshape(p.R_ohm,1,[]); Rb=reshape(p.Rb_ohm,1,[]);
assert(all([L C R Rb p.Cp_F p.Rp_ohm]>0), ...
    'TSCLazo:Pasividad','Se requieren componentes y resistencias positivos.');
assert(p.margenDiDt>=0 && p.margenDiDt<1 && ...
    abs(p.limiteDiDtPlan_A_s-(1-p.margenDiDt)*p.limiteDiDt_A_s) ...
    <=1e-10*max(1,p.limiteDiDt_A_s), ...
    'TSCLazo:Margen','El margen no corresponde al limite de programacion.');

energy=.5*p.Cp_F*a.vCp.^2+sum(.5*L.*a.iRamas.^2+.5*C.*a.vCs.^2,2);
loss=p.Rp_ohm*a.iCp.^2+sum(R.*a.iRamas.^2+a.vCs.^2./Rb,2);
power=a.u.*(a.iCp+sum(a.iRamas,2));
balance=energy-energy(1)-(a.integralEntrada_J-a.integralEntrada_J(1)) ...
    +(a.integralPerdidas_J-a.integralPerdidas_J(1));
balanceTrap=energy-energy(1)-cumtrapz(t,power-loss);
scaleE=max([1;energy]); scaleI=max([1;abs(a.iFuente)]);
scaleV=max([1;abs(a.u);abs(a.vCs(:))]);
q=struct('Muestras',n,'Duracion_s',t(end)-t(1), ...
    'Error_KCL_A',max(abs(a.iFuente-a.iCarga-a.iCp-sum(a.iRamas,2))), ...
    'Error_ESR_Cp_V',max(abs(a.u-a.vCp-p.Rp_ohm*a.iCp)), ...
    'Error_reconstruccion_A',max(abs(a.iCargaReconstruida-a.iCarga)), ...
    'Error_energia_declarada_J',max(abs(energy-a.energia_J)), ...
    'Error_perdidas_declaradas_W',max(abs(loss-a.perdidas_W)), ...
    'Balance_integrales_J',max(abs(balance)), ...
    'Balance_cuadratura_J',max(abs(balanceTrap)), ...
    'I_fuente_RMS_A',sqrt(trapz(t,a.iFuente.^2)/(t(end)-t(1))), ...
    'VC_pico_V',max(abs(a.vCs),[],'all'),'I_rama_pico_A',max(abs(a.iRamas),[],'all'));
q.Tolerancia_balance_integrales_J=1e-5*scaleE;
q.Tolerancia_balance_cuadratura_J=5e-4*max([scaleE,trapz(t,abs(power))+trapz(t,loss)])+1e-6;
q.Corriente_bloqueada_A=max([0;abs(a.iRamas(~logical(a.estadoAplicado)))]);
assert(q.Error_KCL_A<1e-10*scaleI && q.Error_reconstruccion_A<1e-10*scaleI, ...
    'TSCLazo:KCL','No cierra la corriente de nodo o su reconstruccion.');
assert(q.Error_ESR_Cp_V<1e-10*scaleV,'TSCLazo:ESR','La rama Cp no satisface su caida ESR.');
assert(q.Error_energia_declarada_J<1e-10*scaleE && ...
    q.Error_perdidas_declaradas_W<1e-10*max([1;loss]), ...
    'TSCLazo:Magnitudes','Energia o perdidas declaradas inconsistentes.');
assert(q.Balance_integrales_J<q.Tolerancia_balance_integrales_J, ...
    'TSCLazo:Energia','No cierra el balance con los estados integrales.');
assert(q.Balance_cuadratura_J<q.Tolerancia_balance_cuadratura_J, ...
    'TSCLazo:EnergiaCuadratura','No cierra la cuadratura independiente de las trazas.');
assert(q.Corriente_bloqueada_A<1e-8, ...
    'TSCLazo:Bloqueo','Hay corriente no despreciable en una rama bloqueada.');
integrals=auditar_integrales(a,p);
events=auditar_eventos(r);
windows=struct('auditado',false); decisions=struct('auditado',false);
if opciones.AuditarSupervisor
    windows=auditar_ventanas(r); decisions=auditar_decisiones(r);
end
structure=struct('inspeccionado',false,'sin_cronograma_externo',false);
if opciones.InspeccionarModelo, structure=auditar_modelo(r.modelo); end
informe=struct('passed',true,'trazas',q,'leyes_integrales',integrals, ...
    'ventanas',windows,'eventos',events,'decisiones',decisions,'modelo',structure, ...
    'alcance',['Auditoria de la simulacion registrada: causalidad, metrica temporal, ', ...
    'leyes del circuito y maniobras de valvula idealizada. No demuestra superioridad ', ...
    'universal ni extiende el teorema periodico ideal a transitorios.']);
if opciones.Exportar
    folder=opciones.CarpetaSalida;
    if strlength(folder)==0
        folder=fullfile(fileparts(mfilename('fullpath')),'resultados', ...
            'tsc_lazo_cerrado',r.config.CasoTension,r.config.Modo);
    end
    if ~isfolder(folder), mkdir(folder); end
    save(fullfile(folder,'verificacion.mat'),'informe');
    fid=fopen(fullfile(folder,'verificacion.json'),'w'); assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));
    fprintf(fid,'%s\n',jsonencode(informe,PrettyPrint=true));
end
disp(informe);
end

function z=auditar_integrales(a,p)
% Cuadratura local independiente, excluyendo intervalos que contienen una
% maniobra. No diferencia numericamente los estados ni cruza discontinuidades.
t=a.t(:); h=diff(t); on=a.estadoAplicado; ir=a.iRamas; vc=a.vCs;
z=struct('Error_KVL_integral_Vs',zeros(1,2), ...
    'Error_KCL_integral_C',zeros(1,2),'Error_Cp_integral_C',0, ...
    'Intervalos_comprobados',zeros(1,2));
z.Error_Cp_integral_C=max(abs(p.Cp_F*diff(a.vCp)-h.*(a.iCp(1:end-1)+a.iCp(2:end))/2));
for j=1:2
    same=on(1:end-1,j)==on(2:end,j);
    di=on(:,j).*(a.u-p.R_ohm(j)*ir(:,j)-vc(:,j));
    dv=ir(:,j)-vc(:,j)/p.Rb_ohm(j);
    ri=p.L_H(j)*diff(ir(:,j))-h.*(di(1:end-1)+di(2:end))/2;
    rv=p.C_F(j)*diff(vc(:,j))-h.*(dv(1:end-1)+dv(2:end))/2;
    z.Error_KVL_integral_Vs(j)=max([0;abs(ri(same))]);
    z.Error_KCL_integral_C(j)=max([0;abs(rv(same))]);
    z.Intervalos_comprobados(j)=nnz(same);
end
% These thresholds test local equations at the configured solver resolution;
% the report keeps absolute residuals, and refinement quantifies discretization.
z.Tolerancia_KVL_Vs=2e-4*max(1,max(abs(a.u)))*max(h)+1e-8;
z.Tolerancia_KCL_C=2e-4*max(1,max(abs(ir),[],'all'))*max(h)+1e-10;
z.Tolerancia_Cp_C=2e-4*max(1,max(abs(a.iCp)))*max(h)+1e-10;
assert(all(z.Error_KVL_integral_Vs<z.Tolerancia_KVL_Vs), ...
    'TSCLazo:IntegralKVL','La ecuacion integral del inductor no cierra.');
assert(all(z.Error_KCL_integral_C<z.Tolerancia_KCL_C) && ...
    z.Error_Cp_integral_C<z.Tolerancia_Cp_C, ...
    'TSCLazo:IntegralKCL','La ecuacion integral del condensador no cierra.');
end

function z=auditar_ventanas(r)
assert(numel(r.ventanas_decision)==height(r.decisiones), ...
    'TSCLazo:Ventanas','Numero de decisiones y ventanas diferente.');
z=struct('Comprobadas',0,'Anticipo_s',0,'Error_captura_V_A',0, ...
    'Error_reconstruccion_A',0,'Error_coste_RMS_A2',0,'Error_Gram_A2',0, ...
    'Error_coste_control_A2',0,'Residuo_KVL_V',0,'Residuo_KCL_A',0, ...
    'Residuo_potencia_W',0);
p=r.parametros; tr=r.trazas; d=r.decisiones; bits=p.bits;
for k=1:numel(r.ventanas_decision)
    a=r.ventanas_decision{k}; if isempty(a), continue; end
    z.Comprobadas=z.Comprobadas+1; N=numel(a.t); tw=a.t(:);
    z.Anticipo_s=max([z.Anticipo_s;tw-d.Tiempo_s(k)]);
    assert(all(diff(tw)>0) && max(tw)<d.Tiempo_s(k)-1e-12 && ...
        abs(tw(1)-(d.Tiempo_s(k)-p.T))<1e-10, ...
        'TSCLazo:VentanaCausal','La ventana no es el ultimo periodo adquirido completo.');
    assert(max(abs(diff(tw)-p.T/N))<1e-11, ...
        'TSCLazo:Muestreo','La ventana no tiene paso temporal uniforme.');
    assert(min(tw)>=tr.t(1)-1e-10 && max(tw)<=tr.t(end)+1e-10, ...
        'TSCLazo:Adquisicion','La ventana queda fuera del intervalo registrado.');
    % Recorta solo el redondeo de un extremo temporal ya validado, nunca
    % extrapola una medida que aun no haya sido adquirida.
    query=min(max(tw,tr.t(1)),tr.t(end));
    idx=interp1(tr.t,(1:numel(tr.t))',query,'nearest');
    assert(all(isfinite(idx)) && max(abs(tr.t(idx)-tw))<1e-10, ...
        'TSCLazo:Adquisicion','Faltan las muestras originales en las trazas.');
    measured=[a.u,a.iFuente,a.iCp,a.iRamas];
    source=[tr.u(idx),tr.iFuente(idx),tr.iCp(idx),tr.iRamas(idx,:)];
    z.Error_captura_V_A=max(z.Error_captura_V_A,max(abs(measured-source),[],'all'));
    rec=a.iFuente-a.iCp-sum(a.iRamas,2);
    z.Error_reconstruccion_A=max([z.Error_reconstruccion_A; ...
        abs(rec-a.iCargaReconstruida);abs(rec-tr.iCarga(idx))]);
    base=rec+a.iCpReferencia; phi=a.phiRamas;
    direct=mean((base+phi*bits').^2,1);
    z.Error_coste_RMS_A2=max(z.Error_coste_RMS_A2,max(abs(direct-a.costesRMS(:)')));
    gram=phi'*phi/N; lin=phi'*base/N; constant=mean(base.^2);
    z.Error_Gram_A2=max([z.Error_Gram_A2;abs(gram(:)-a.G(:)); ...
        abs(lin(:)-a.h(:));abs(constant-a.c)]);
    expected=direct;
    if r.config.Modo=="convencional"
        co=cos(p.omega*tw); sn=sin(p.omega*tw);
        ac=2*mean(a.u.*co); as=2*mean(a.u.*sn);
        quadrature=sqrt(2)*(ac*sn-as*co)/hypot(ac,as);
        expected=mean((base+phi*bits').*quadrature,1).^2;
    end
    z.Error_coste_control_A2=max(z.Error_coste_control_A2,max(abs(expected-a.costesControl(:)')));
    Dt=p.omega*derivador_periodico(N); vcp=a.u-p.Rp_ohm*a.iCpReferencia;
    z.Residuo_KCL_A=max(z.Residuo_KCL_A,max(abs(a.iCpReferencia-p.Cp_F*Dt*vcp)));
    z.Residuo_potencia_W=max(z.Residuo_potencia_W, ...
        abs(mean(a.u.*a.iCpReferencia)-p.Rp_ohm*mean(a.iCpReferencia.^2)));
    for j=1:2
        state=a.referenciaEstado{j}; ii=state(:,1); vv=state(:,2);
        assert(max(abs(ii-phi(:,j)))<1e-10,'TSCLazo:PlantillaEstado','La plantilla y su estado difieren.');
        z.Residuo_KVL_V=max(z.Residuo_KVL_V,max(abs(a.u-p.R_ohm(j)*ii-vv-p.L_H(j)*Dt*ii)));
        z.Residuo_KCL_A=max(z.Residuo_KCL_A,max(abs(ii-vv/p.Rb_ohm(j)-p.C_F(j)*Dt*vv)));
        z.Residuo_potencia_W=max(z.Residuo_potencia_W, ...
            abs(mean(a.u.*ii)-p.R_ohm(j)*mean(ii.^2)-mean(vv.^2)/p.Rb_ohm(j)));
    end
end
assert(z.Comprobadas>0,'TSCLazo:SinDecision','No hay decisiones auditables.');
assert(z.Anticipo_s<1e-10 && z.Error_captura_V_A<1e-6 && z.Error_reconstruccion_A<1e-7, ...
    'TSCLazo:Causalidad','La ventana no corresponde a medidas efectivamente adquiridas.');
assert(max([z.Error_coste_RMS_A2,z.Error_Gram_A2,z.Error_coste_control_A2])<1e-8, ...
    'TSCLazo:Gram','El coste no coincide con el producto temporal directo.');
assert(z.Residuo_KVL_V<1e-7 && z.Residuo_KCL_A<1e-8 && z.Residuo_potencia_W<1e-7, ...
    'TSCLazo:Referencias','Las referencias periodicas no satisfacen sus leyes de circuito.');
end

function z=auditar_eventos(r)
e=r.eventos; a=r.trazas; p=r.parametros;
z=struct('Maniobras',0,'Rechazos',0,'Residuo_OFF_A',0,'Energia_reset_J',0, ...
    'Error_VC_evento_V',0,'Max_di_dt_ON_A_s',0,'Error_di_dt_ON_A_s',0, ...
    'Min_retardo_orden_s',Inf,'Error_instante_programado_s',0);
for k=1:height(e)
    j=e.Rama(k); ta=e.Tiempo_s(k); kind=string(e.Tipo(k));
    assert(ismember(j,[1 2]),'TSCLazo:EventoRama','Rama de evento incorrecta.');
    [distance,ix]=min(abs(a.t-ta));
    assert(distance<1e-8, ...
        'TSCLazo:EventoTiempo','No hay traza de planta en el instante de maniobra.');
    z.Error_VC_evento_V=max(z.Error_VC_evento_V,abs(a.vCs(ix,j)-e.VC_pre_V(k)));
    if kind=="rechazo_ON", z.Rechazos=z.Rechazos+1; continue; end
    assert(ismember(kind,["ON","OFF"]),'TSCLazo:TipoEvento','Tipo de evento desconocido.');
    z.Maniobras=z.Maniobras+1;
    target=double(kind=="ON"); requested=a.estadoSolicitado(1:ix,j);
    changes=find(diff([0;requested])~=0 & requested==target);
    assert(~isempty(changes),'TSCLazo:EventoSinOrden','Una maniobra no tiene orden antecedente.');
    order=a.t(changes(end)); z.Min_retardo_orden_s=min(z.Min_retardo_orden_s,ta-order);
    assert(order<=ta+1e-10 && requested(end)==target, ...
        'TSCLazo:OrdenCausal','La maniobra no respeta la orden vigente.');
    assert(a.estadoAplicado(ix,j)==target,'TSCLazo:Aplicacion','Evento y estado aplicado discrepan.');
    if kind=="OFF"
        z.Residuo_OFF_A=max(z.Residuo_OFF_A,abs(e.I_pre_A(k)));
        z.Energia_reset_J=z.Energia_reset_J+.5*p.L_H(j)*e.I_pre_A(k)^2;
    else
        slope=(a.u(ix)-e.VC_pre_V(k))/p.L_H(j);
        z.Max_di_dt_ON_A_s=max(z.Max_di_dt_ON_A_s,abs(slope));
        z.Error_di_dt_ON_A_s=max(z.Error_di_dt_ON_A_s,abs(slope-e.DiDt_A_s(k)));
        assert(e.Admisible(k) && abs(e.I_pre_A(k))<1e-8, ...
            'TSCLazo:Conexion','Se conecta una rama no admisible o con corriente inicial.');
        plans=r.decisiones.(sprintf('ProgramacionON%d_s',j));
        [delay,idxPlan]=min(abs(plans-ta));
        assert(delay<1e-8 && r.decisiones.Tiempo_s(idxPlan)<=ta+1e-10, ...
            'TSCLazo:ProgramaCausal','El ON no corresponde a una programacion antecedente.');
        z.Error_instante_programado_s=max(z.Error_instante_programado_s,delay);
        predicted=r.decisiones.(sprintf('DiDtPredicha%d_A_s',j));
        assert(isfinite(predicted(idxPlan)) && abs(predicted(idxPlan)) ...
            <=p.limiteDiDtPlan_A_s+1e-7*max(1,p.limiteDiDtPlan_A_s), ...
            'TSCLazo:MargenPlan','La programacion no respeta el margen de pendiente.');
    end
end
for j=1:2
    changes=diff([0;a.estadoAplicado(:,j)]);
    assert(nnz(changes>0)==nnz(e.Rama==j & e.Tipo=="ON") && ...
        nnz(changes<0)==nnz(e.Rama==j & e.Tipo=="OFF"), ...
        'TSCLazo:RegistroManiobras','Los cambios de conduccion y eventos no coinciden.');
end
z.Tolerancia_OFF_A=1e-4; z.Tolerancia_energia_reset_J=1e-8;
assert(z.Residuo_OFF_A<z.Tolerancia_OFF_A && z.Energia_reset_J<z.Tolerancia_energia_reset_J, ...
    'TSCLazo:Extincion','La apertura elimina corriente o energia no despreciable.');
assert(z.Error_VC_evento_V<1e-5, ...
    'TSCLazo:ContinuidadVC','VC registrado antes de maniobrar difiere de la planta posterior.');
assert(z.Max_di_dt_ON_A_s<=p.limiteDiDt_A_s+1e-7*max(1,p.limiteDiDt_A_s) && ...
    z.Error_di_dt_ON_A_s<1e-5, ...
    'TSCLazo:PendienteON','La pendiente real incumple el limite o el registro.');
end

function z=auditar_decisiones(r)
d=r.decisiones; e=r.eventos; p=r.parametros;
z=struct('Numero',height(d),'Anticipo_s',0,'Permanencia_minima_s',Inf, ...
    'Cambios',0,'Error_ganancia_A2',0);
if isempty(d), return; end
z.Anticipo_s=max([0;d.MaxMuestraUsada_s-d.Tiempo_s]);
assert(z.Anticipo_s<1e-10,'TSCLazo:DecisionFutura','La decision declara muestras futuras.');
assert(all(~logical(d.Cambio)|~logical(d.Pendiente)), ...
    'TSCLazo:OrdenPendiente','Se sustituye una orden aun pendiente.');
for k=1:height(d)
    costs=r.ventanas_decision{k}.costesControl(:);
    actual=d.EstadoActual(k)+1; candidate=d.EstadoCandidato(k)+1;
    assert(ismember(actual,1:4) && ismember(candidate,1:4), ...
        'TSCLazo:Candidato','Indice de combinacion fuera del diccionario.');
    if r.config.Modo=="fijo"
        expected=find(all(p.bits==r.config.BitsFijos,2));
    else
        best=min(costs); ties=find(abs(costs-best)<1e-10*max(1,best));
        [~,jTie]=min(sum(abs(p.bits(ties,:)-p.bits(actual,:)),2));
        expected=ties(jTie);
    end
    assert(candidate==expected,'TSCLazo:Seleccion','No se elige el candidato declarado por el criterio.');
    gain=costs(actual)-costs(candidate);
    z.Error_ganancia_A2=max(z.Error_ganancia_A2,abs(gain-d.Ganancia_A2(k)));
    allowed=~d.Pendiente(k) && d.Permanencia_s(k)>=p.permanenciaMin_s-1e-12 && ...
        (gain>p.histeresis_A2 || r.config.Modo=="fijo");
    assert(logical(d.Cambio(k))==(allowed && candidate~=actual), ...
        'TSCLazo:Histeresis','La decision incumple histeresis, permanencia o bloqueo pendiente.');
end
assert(z.Error_ganancia_A2<1e-8,'TSCLazo:Ganancia','Ganancia declarada incorrecta.');
for k=find(logical(d.Cambio))'
    z.Cambios=z.Cambios+1;
    previous=e.Tiempo_s(ismember(e.Tipo,["ON","OFF"]) & e.Tiempo_s<d.Tiempo_s(k)-1e-10);
    if isempty(previous), continue; end
    age=d.Tiempo_s(k)-max(previous);
    z.Permanencia_minima_s=min(z.Permanencia_minima_s,age);
    assert(age>=p.permanenciaMin_s-1e-9, ...
        'TSCLazo:Permanencia','Se incumple la permanencia desde la ultima maniobra real.');
end
end

function z=auditar_modelo(path)
assert(isfile(path),'TSCLazo:Modelo','No existe el modelo guardado.');
[~,name]=fileparts(path); wasLoaded=bdIsLoaded(name);
if ~wasLoaded, load_system(path); end
cleanup=onCleanup(@()cerrar_si_procede(name,wasLoaded));
sources=find_system(name,'LookUnderMasks','all','FollowLinks','on','BlockType','FromWorkspace');
variables=cellfun(@(b)get_param(b,'VariableName'),sources,'UniformOutput',false);
assert(all(strcmp(variables,'tm_carga')), ...
    'TSCLazo:CronogramaExterno','Se detecta una entrada externa distinta de la carga exogena.');
supervisor=[name,'/Supervisor_temporal'];
assert(strcmp(get_param(supervisor,'FunctionName'),'tsc_supervisor_sfun'), ...
    'TSCLazo:Supervisor','No se encuentra el supervisor ejecutado dentro del modelo.');
assert(strcmp(get_param([name,'/Feedback_causal'],'BlockType'),'Memory'), ...
    'TSCLazo:Feedback','El lazo debe declarar la memoria usada para evitar feedthrough.');
for j=1:2
    branch=[name,sprintf('/LC_%d',j)];
    assert(strcmp(get_param([branch,'/Estado_VC'],'ExternalReset'),'none'), ...
        'TSCLazo:ResetVC','El capacitor tiene un reset externo.');
    assert(strcmp(get_param([branch,'/Estado_I'],'ExternalReset'),'falling'), ...
        'TSCLazo:ResetI','El reset de corriente no se limita al paso de conduccion a bloqueo.');
    assert(strcmp(get_param([branch,'/Estado_I'],'ShowStatePort'),'on'), ...
        'TSCLazo:CorrientePrevia','Falta la salida de estado anterior al reset del inductor.');
    executor=[name,sprintf('/Ejecutor_%d',j)];
    assert(strcmp(get_param([executor,'/Latch_evento/Actualizar'],'FunctionName'), ...
        'tsc_ejecutor_sfun'),'TSCLazo:Ejecutor','Falta el ejecutor interno de maniobras.');
    connections=get_param(executor,'PortConnectivity');
    assert(connections(5).SrcPort==2 && ...
        strcmp(get_param(connections(5).SrcBlock,'Name'),sprintf('LC_%d',j)), ...
        'TSCLazo:CorrientePrevia','El ejecutor no recibe el estado previo al reset.');
    assert(strcmp(get_param([name,sprintf('/Ejecutor_%d/Cero_corriente',j)], ...
        'BlockType'),'HitCross'), ...
        'TSCLazo:CeroContinuo','No se detecta el cero de corriente en el dominio continuo.');
end
ports=get_param(supervisor,'PortConnectivity');
assert(numel(ports)>=3 && ports(1).SrcBlock>0 && ports(2).SrcBlock>0 && ports(3).SrcBlock>0, ...
    'TSCLazo:SupervisorEntradas','Faltan medidas o feedback del supervisor.');
assert(strcmp(get_param(ports(1).SrcBlock,'Name'),'Medidas') && ...
    strcmp(get_param(ports(2).SrcBlock,'Name'),'Feedback_causal') && ...
    strcmp(get_param(ports(3).SrcBlock,'Name'),'VC_actuales'), ...
    'TSCLazo:SupervisorConectado','El supervisor no lee medidas, feedback y VC de la planta.');
z=struct('inspeccionado',true,'sin_cronograma_externo',true, ...
    'entradas_workspace',{variables},'supervisor','tsc_supervisor_sfun', ...
    'feedback','Memory: un paso mayor del solver', ...
    'adquisicion_y_VC_separadas',true, ...
    'VC_sin_reset',true,'OFF_deteccion_continua',true,'I_pre_desde_estado_inductor',true, ...
    'alcance','Inspeccion estructural del SLX guardado; se complementa con trazas de ejecucion.');
end

function cerrar_si_procede(name,wasLoaded)
if ~wasLoaded && bdIsLoaded(name), close_system(name,0); end
end
