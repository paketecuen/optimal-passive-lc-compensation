function informe = verificar_control_tsc(datos, opciones)
%VERIFICAR_CONTROL_TSC Auditoria de trazas, eventos y decisiones del banco TSC.
% No vuelve a simular: recibe un registro, una campana o carga el MAT guardado.
% El balance por cuadratura se calcula de las magnitudes fisicas registradas,
% independientemente de los estados integrales usados por el motor ODE.
arguments
    datos = []
    opciones.Exportar (1,1) logical = true
    opciones.Refinado = []
end
carpeta=fullfile(fileparts(mfilename('fullpath')),'resultados','control_tsc');
if isempty(datos)
    s=load(fullfile(carpeta,'control_tsc.mat'),'saved'); datos=s.saved;
end
if isfield(datos,'registros'), registros=datos.registros; else, registros={datos}; end
filas=struct([]);
for k=1:numel(registros)
    r=registros{k}; a=r.trazas; p=r.parametros; e=r.eventos; d=r.decisiones;
    assert(all(isfield(p,{'limiteDiDt_A_s','permanenciaMin_s','T'})), ...
        'TSC:Contrato','Faltan limites o periodo en los parametros auditables.');
    assert(isfield(a,'iCargaReconstruida') && isfield(r,'ventanas_decision'), ...
        'TSC:Contrato','Faltan reconstruccion o ventanas para una auditoria independiente.');
    limitePlan=p.limiteDiDt_A_s; margen=0;
    if isfield(p,'margenDiDt')
        assert(isfield(p,'limiteDiDtPlan_A_s'),'TSC:Margen','Falta el limite de programacion.');
        margen=p.margenDiDt; limitePlan=p.limiteDiDtPlan_A_s;
        assert(margen>=0 && margen<1 && limitePlan<=p.limiteDiDt_A_s, ...
            'TSC:Margen','Margen de prediccion incorrecto.');
        if isfinite(p.limiteDiDt_A_s)
            assert(abs(limitePlan-(1-margen)*p.limiteDiDt_A_s)<1e-10*max(1,p.limiteDiDt_A_s), ...
                'TSC:Margen','El margen altera el limite fisico en vez del limite de programacion.');
        end
    end
    t=a.t(:); n=numel(t); assert(n>2 && all(diff(t)>0),'TSC:Tiempo','Tiempo no creciente.');
    assert(size(a.iRamas,2)==2 && size(a.vCs,2)==2,'TSC:Ramas','Se esperan dos LC conmutables.');
    assert(all(isfinite([a.u(:);a.iFuente(:);a.iCarga(:);a.vCs(:);a.iRamas(:)])), ...
        'TSC:Finito','Registro no finito.');
    L=reshape(p.L_H,1,[]); C=reshape(p.C_F,1,[]);
    R=reshape(p.R_ohm,1,[]); Rb=reshape(p.Rb_ohm,1,[]);
    energia=.5*p.Cp_F*a.vCp(:).^2+sum(.5*L.*a.iRamas.^2+.5*C.*a.vCs.^2,2);
    perdidas=p.Rp_ohm*a.iCp(:).^2+sum(R.*a.iRamas.^2+a.vCs.^2./Rb,2);
    entrada=a.u(:).*(a.iCp(:)+sum(a.iRamas,2));
    energiaError=max(abs(energia-a.energia_J(:)));
    perdidaError=max(abs(perdidas-a.perdidas_W(:)));
    balanceODE=max(abs(energia-energia(1)-a.integralEntrada_J(:)+a.integralPerdidas_J(:)));
    balanceCuadratura=max(abs(energia-energia(1)-cumtrapz(t,entrada-perdidas)));
    escalaEnergia=max([max(energia),trapz(t,abs(entrada))+trapz(t,perdidas),1e-6]);
    toleranciaCuadratura=5e-4*escalaEnergia+1e-6;
    assert(energiaError<1e-9*max(1,max(energia)),'TSC:Energia','Energia declarada incorrecta.');
    assert(perdidaError<1e-9*max(1,max(perdidas)),'TSC:Perdidas','Perdidas declaradas incorrectas.');
    assert(balanceODE<1e-5*max(1,max(energia)),'TSC:BalanceODE','El balance del motor no cierra.');
    assert(balanceCuadratura<toleranciaCuadratura,'TSC:BalanceIndependiente', ...
        'No cierra el balance mediante cuadratura independiente sobre las trazas.');
    kirchhoff=max(abs(a.iFuente(:)-a.iCarga(:)-a.iCp(:)-sum(a.iRamas,2)));
    assert(kirchhoff<1e-10*max(1,max(abs(a.iFuente))),'TSC:Kirchhoff','No cierra la corriente de nodo.');
    assert(all(ismember(a.estadoAplicado(:),[0,1])) && all(ismember(a.estadoSolicitado(:),[0,1])), ...
        'TSC:Estados','Estados binarios incorrectos.');
    assert(all(ismember(a.pendiente(:),[-1,0,1])),'TSC:Pendiente','Estado pendiente incorrecto.');
    assert(~any(a.pendiente(:)==1 & logical(a.estadoAplicado(:))) && ...
        ~any(a.pendiente(:)==-1 & ~logical(a.estadoAplicado(:))), ...
        'TSC:PendienteAplicado','El pendiente no concuerda con la conduccion actual.');
    bloqueada=~logical(a.estadoAplicado);
    corrienteBloqueada=max(abs(a.iRamas(bloqueada)),[],'all');
    if isempty(corrienteBloqueada), corrienteBloqueada=0; end
    assert(corrienteBloqueada<1e-8,'TSC:Bloqueo','Una rama bloqueada conserva corriente.');
    [eventos,confirmados]=auditar_eventos(e,p);
    decisiones=auditar_decisiones(d,e,confirmados,p);
    reconstruccion=auditar_reconstruccion(a);
    gram=auditar_gram(r);
    filas=[filas;struct('Registro',k,'Muestras',n,'Eventos',height(e), ...
        'Decisiones',height(d),'Error_Kirchhoff_A',kirchhoff, ...
        'Error_energia_J',energiaError,'Error_perdidas_W',perdidaError, ...
        'Balance_ODE_J',balanceODE,'Balance_cuadratura_J',balanceCuadratura, ...
        'Tolerancia_cuadratura_J',toleranciaCuadratura, ...
        'Corriente_bloqueada_A',corrienteBloqueada, ...
        'Salto_VC_V',eventos.saltoVC,'Residuo_OFF_A',eventos.residuoOFF, ...
        'Max_di_dt_ON_A_s',eventos.maxDiDt,'Limite_di_dt_real_A_s',p.limiteDiDt_A_s, ...
        'Limite_di_dt_plan_A_s',limitePlan,'Margen_prediccion',margen, ...
        'Rechazos_ON',eventos.rechazos, ...
        'Timeouts_OFF',eventos.timeouts,'Causalidad_s',decisiones.anticipo, ...
        'Permanencia_minima_s',decisiones.minDwell, ...
        'Error_reconstruccion_A',reconstruccion,'Error_Gram_A2',gram.coste, ...
        'Residuo_plantillas_KVL_V',gram.kvl,'Residuo_plantillas_KCL_A',gram.kcl, ...
        'Balance_plantillas_W',gram.potencia)]; %#ok<AGROW>
end
tabla=struct2table(filas);
informe=struct('passed',true,'registros',numel(registros),'resultados',tabla, ...
    'alcance','Auditoria numerica del modelo hibrido; no valida un SCR ni hardware.');
if ~isempty(opciones.Refinado)
    fino=verificar_control_tsc(opciones.Refinado,Exportar=false);
    informe.convergencia=comparar_refinado(registros{1},opciones.Refinado,tabla(1,:),fino.resultados);
end
if opciones.Exportar
    if ~isfolder(carpeta), mkdir(carpeta); end
    writetable(tabla,fullfile(carpeta,'verificacion.csv'));
    copia=informe; copia.resultados=table2struct(tabla);
    fid=fopen(fullfile(carpeta,'verificacion.json'),'w'); assert(fid>=0);
    c=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(copia,PrettyPrint=true));
end
disp(tabla);
end

function c=comparar_refinado(base,fino,verBase,verFino)
assert(base.parametros.pasoMuestreo_s>fino.parametros.pasoMuestreo_s, ...
    'TSC:Refinamiento','El segundo registro debe tener paso de salida menor.');
opcionesBase=rmfield(base.parametros.opciones,'PasoSalida_s');
opcionesFino=rmfield(fino.parametros.opciones,'PasoSalida_s');
assert(isequaln(orderfields(opcionesBase),orderfields(opcionesFino)),'TSC:Refinamiento', ...
    'Las estrategias, carga, condiciones iniciales y limites deben coincidir.');
assert(base.parametros.PNominal_W==fino.parametros.PNominal_W && ...
    isequal(base.parametros.coefTension,fino.parametros.coefTension) && ...
    base.parametros.limiteDiDt_A_s==fino.parametros.limiteDiDt_A_s && ...
    base.parametros.margenDiDt==fino.parametros.margenDiDt && ...
    abs(base.trazas.t(end)-fino.trazas.t(end))<1e-10, ...
    'TSC:Refinamiento','Se comparan escenarios o limites distintos.');
a=base.eventos(ismember(base.eventos.Tipo,["ON","OFF"]),:);
b=fino.eventos(ismember(fino.eventos.Tipo,["ON","OFF"]),:);
mismaSecuencia=height(a)==height(b) && isequal(a.Tipo,b.Tipo) && isequal(a.Rama,b.Rama);
deltaEventos=NaN;
if mismaSecuencia
    if isempty(a), deltaEventos=0; else, deltaEventos=max(abs(a.Tiempo_s-b.Tiempo_s)); end
end
mismaCronologia=height(base.decisiones)==height(fino.decisiones) && ...
    max(abs(base.decisiones.Tiempo_s-fino.decisiones.Tiempo_s))<1e-10;
decisionesDistintas=NaN;
if mismaCronologia
    decisionesDistintas=sum(base.decisiones.EstadoSolicitado~=fino.decisiones.EstadoSolicitado);
end
irms0=sqrt(trapz(base.trazas.t,base.trazas.iFuente.^2)/base.trazas.t(end));
irms1=sqrt(trapz(fino.trazas.t,fino.trazas.iFuente.^2)/fino.trazas.t(end));
c=struct('paso_base_s',base.parametros.pasoMuestreo_s, ...
    'paso_refinado_s',fino.parametros.pasoMuestreo_s, ...
    'I_RMS_base_A',irms0,'I_RMS_refinado_A',irms1, ...
    'diferencia_I_RMS_A',abs(irms1-irms0), ...
    'diferencia_I_RMS_relativa',abs(irms1-irms0)/max(irms1,eps), ...
    'diferencia_pico_corriente_A',abs(max(abs(base.trazas.iRamas),[],'all')-max(abs(fino.trazas.iRamas),[],'all')), ...
    'diferencia_pico_VC_V',abs(max(abs(base.trazas.vCs),[],'all')-max(abs(fino.trazas.vCs),[],'all')), ...
    'secuencia_maniobras_igual',mismaSecuencia,'diferencia_instantes_maniobra_s',deltaEventos, ...
    'cronologia_decisiones_igual',mismaCronologia,'decisiones_distintas',decisionesDistintas, ...
    'rechazos_base',sum(base.eventos.Tipo=="rechazo_ON"), ...
    'rechazos_refinado',sum(fino.eventos.Tipo=="rechazo_ON"), ...
    'balance_cuadratura_base_J',verBase.Balance_cuadratura_J(1), ...
    'balance_cuadratura_refinado_J',verFino.Balance_cuadratura_J(1), ...
    'interpretacion','Diferencias cuantificadas; igualdad de maniobras no sustituye validacion de hardware.');
disp(c);
end

function [a,confirmados]=auditar_eventos(e,p)
a=struct('saltoVC',0,'residuoOFF',0,'maxDiDt',0,'rechazos',0,'timeouts',0);
confirmados=false(height(e),1);
if isempty(e), return; end
tipo=lower(string(e.Tipo));
on=ismember(tipo,["on","conectar","conexion"]);
off=ismember(tipo,["off","bloquear","desconexion"]);
confirmados=on|off;
a.rechazos=sum(tipo=="rechazo_on"); a.timeouts=sum(tipo=="timeout_off");
assert(all(e.Tiempo_s>=e.Orden_s-1e-10),'TSC:EventoCausal','Maniobra anterior a orden.');
assert(max(abs(e.Retardo_s-(e.Tiempo_s-e.Orden_s)))<1e-9,'TSC:Retardo','Retardo inconsistente.');
if any(confirmados)
    a.saltoVC=max(abs(e.VC_post_V(confirmados)-e.VC_pre_V(confirmados)));
    assert(a.saltoVC<1e-9,'TSC:Continuidad','Se reinicia la tension del condensador.');
end
if any(off)
    a.residuoOFF=max(abs([e.I_pre_A(off);e.I_post_A(off)]));
    assert(a.residuoOFF<1e-7,'TSC:Extincion','Se abre una rama fuera del cero de corriente.');
end
if any(on)
    assert(all(e.Admisible(on)),'TSC:Admisible','Se confirma una conexion no admisible.');
    assert(max(abs([e.I_pre_A(on);e.I_post_A(on)]))<1e-8,'TSC:ONCorriente','Conexion con corriente inicial no nula.');
    a.maxDiDt=max(abs(e.DiDt_A_s(on)));
    limite=p.limiteDiDt_A_s;
    assert(a.maxDiDt<=limite+1e-7*max(1,limite)+1e-8,'TSC:Limite','Se supera el limite de pendiente al conectar.');
    if isfield(p,'margenDiDt')
        assert(ismember('DiDt_predicha_A_s',e.Properties.VariableNames), ...
            'TSC:Prediccion','Falta la pendiente usada para programar la maniobra.');
        predicha=e.DiDt_predicha_A_s(on); plan=p.limiteDiDtPlan_A_s;
        assert(all(isfinite(predicha)) && all(abs(predicha)<=plan+1e-7*max(1,plan)+1e-8), ...
            'TSC:Prediccion','La maniobra programada no respeta su margen.');
    end
end
end

function a=auditar_decisiones(d,e,confirmados,p)
a=struct('anticipo',0,'minDwell',Inf);
if isempty(d), return; end
a.anticipo=max([0;d.FinVentana_s-d.Tiempo_s;d.MaxMuestraUsada_s-d.Tiempo_s]);
assert(a.anticipo<1e-10,'TSC:Causalidad','Se emplean muestras futuras.');
assert(all(d.InicioVentana_s>=-1e-10) && all(d.InicioVentana_s<d.FinVentana_s), ...
    'TSC:Ventana','Ventana temporal incorrecta.');
% Interpolar una muestra interior puede necesitar el vecino adquirido mas
% tarde que FinVentana_s; sigue siendo causal si no supera Tiempo_s.
assert(all(~logical(d.Cambio)|~logical(d.Pendiente)),'TSC:OrdenPendiente', ...
    'Se sustituye una orden mientras otra esta pendiente.');
limite=p.permanenciaMin_s;
for j=find(logical(d.Cambio))'
    % Una maniobra inmediata nacida de esta misma decision ocurre despues
    % de evaluarla; no debe convertirse retrospectivamente en su permanencia.
    anteriores=e.Tiempo_s(confirmados & e.Tiempo_s<=d.Tiempo_s(j)+1e-10 & ...
        e.Orden_s<d.Tiempo_s(j)-1e-10);
    if isempty(anteriores), continue; end
    intervalo=d.Tiempo_s(j)-max(anteriores);
    a.minDwell=min(a.minDwell,intervalo);
    assert(intervalo>=limite-1e-9,'TSC:Permanencia', ...
        'La permanencia debe contarse desde la ultima maniobra confirmada.');
end
end

function error=auditar_reconstruccion(a)
error=NaN;
if isfield(a,'iCargaReconstruida')
    cp=a.iCp; ramas=a.iRamas;
    if isfield(a,'iCpMedida'), cp=a.iCpMedida; end
    if isfield(a,'iRamasMedidas'), ramas=a.iRamasMedidas; end
    fuente=a.iFuente;
    if isfield(a,'iFuenteMedida'), fuente=a.iFuenteMedida; end
    rec=fuente(:)-cp(:)-sum(ramas,2);
    error=max(abs(rec-a.iCargaReconstruida(:)));
    assert(error<1e-10*max(1,max(abs(rec))),'TSC:Reconstruccion', ...
        'La carga no se reconstruye con las corrientes realmente medidas.');
end
end

function error=auditar_gram(r)
assert(~isempty(r.ventanas_decision),'TSC:Ventanas','No hay ventanas para auditar el coste.');
error=struct('coste',0,'kvl',0,'kcl',0,'potencia',0); comprobadas=0;
p=r.parametros;
for j=1:numel(r.ventanas_decision)
    a=r.ventanas_decision{j};
    if isempty(a), continue; end
    comprobadas=comprobadas+1;
    assert(isfield(a,'costesRMS'),'TSC:Contrato','Faltan costes RMS declarados.');
    rec=a.iFuente(:)-a.iCp(:)-sum(a.iRamas,2);
    assert(max(abs(rec-a.iCargaReconstruida(:)))<1e-10*max(1,max(abs(rec))), ...
        'TSC:ReconstruccionVentana','La ventana resta plantillas en vez de corrientes reales.');
    assert(max(a.t)<=r.decisiones.Tiempo_s(j)+1e-10, ...
        'TSC:VentanaFutura','La ventana guardada contiene muestras futuras.');
    bits=[0,0;0,1;1,0;1,1];
    if isfield(r.parametros,'bits'), bits=r.parametros.bits; end
    base=a.iCargaReconstruida(:)+a.iCpReferencia(:);
    directo=mean((base+a.phiRamas*bits').^2,1);
    costes=a.costesRMS(:)';
    error.coste=max(error.coste,max(abs(directo-costes)));
    % Residuales de las leyes electricas en SI, sin volver a resolver el
    % sistema usado para obtener las plantillas del supervisor.
    Dt=p.omega*derivador_periodico(numel(a.t));
    vCp=a.u-p.Rp_ohm*a.iCpReferencia;
    error.kcl=max(error.kcl,max(abs(a.iCpReferencia-p.Cp_F*Dt*vCp)));
    error.potencia=max(error.potencia,abs(mean(a.u.*a.iCpReferencia)- ...
        p.Rp_ohm*mean(a.iCpReferencia.^2)));
    for rama=1:2
        estado=a.referenciaEstado{rama}; i=estado(:,1); v=estado(:,2);
        error.kvl=max(error.kvl,max(abs(a.u-p.R_ohm(rama)*i-v-p.L_H(rama)*Dt*i)));
        error.kcl=max(error.kcl,max(abs(i-v/p.Rb_ohm(rama)-p.C_F(rama)*Dt*v)));
        error.potencia=max(error.potencia,abs(mean(a.u.*i)- ...
            p.R_ohm(rama)*mean(i.^2)-mean(v.^2)/p.Rb_ohm(rama)));
    end
end
assert(comprobadas>0,'TSC:Ventanas','Las ventanas guardadas estan vacias.');
assert(error.coste<1e-8,'TSC:Gram','El coste geometrico difiere de la suma temporal directa.');
assert(error.kvl<1e-7 && error.kcl<1e-8 && error.potencia<1e-7, ...
    'TSC:Plantillas','Las plantillas no satisfacen las leyes de los componentes con perdidas.');
end
