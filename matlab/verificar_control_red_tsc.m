function report = verificar_control_red_tsc(r,opciones)
%VERIFICAR_CONTROL_RED_TSC Circuito, prediccion, admision y causalidad registradas.
% No llama al supervisor ni al planificador para recomputar sus decisiones.
arguments
    r (1,1) struct
    opciones.InspeccionarModelo (1,1) logical = true
end
base=verificar_fuente_tsc(r,InspeccionarModelo=false,AuditarSupervisor=false);
cfg=r.config; p=cfg.p; c=cfg.ControlRed; a=r.trazas; d=r.decisiones;
D=p.omega*derivador_periodico(63); bits=p.bits;
q=struct('Decisiones',height(d),'Error_captura',0,'Residuo_periodico',0, ...
    'Error_costes_A2',0,'Error_estado_actual',0,'ON_admitidos',0,'ON_rechazados',0, ...
    'Error_cota_evento',0,'Error_energia_evento_J',0,'Rechazos_plan',sum(d.RechazosPlan));
assert(height(d)>0,'TSCRed:Decisiones','No hay decisiones auditables.');
for k=1:height(d)
    w=r.ventanas_decision{k}; tw=w.t; td=d.Tiempo_s(k);
    assert(max(tw)<td-1e-12 && abs(tw(1)-td+p.T)<1e-10 && ...
        max(abs(diff(tw)-cfg.TsSupervision))<1e-10, ...
        'TSCRed:Causalidad','Ventana no causal o muestreo incorrecto.');
    ix=interp1(a.t,(1:numel(a.t))',max(0,tw),'nearest');
    assert(max(abs(a.t(ix)-tw))<1e-9,'TSCRed:Captura','Faltan muestras originales.');
    expected=[a.u(ix),a.iFuente(ix),a.iCp(ix),a.iRamas(ix,:),a.vCs(ix,:),a.e(ix)];
    q.Error_captura=max(q.Error_captura,max(abs(w.medidas-expected),[],'all'));
    il=w.medidas(:,2)-w.medidas(:,3)-sum(w.medidas(:,4:5),2);
    assert(max(abs(il-w.iCargaReconstruida))<1e-10,'TSCRed:Reconstruccion','Carga incorrecta.');
    costs=zeros(1,4); frozen=costs;
    I=eye(63); icpFixed=(I+p.Rp_ohm*p.Cp_F*D)\(p.Cp_F*D*w.u);
    phi=zeros(63,2);
    for j=1:2
        h=[p.L_H(j)*D+p.R_ohm(j)*I,I;-I,p.C_F(j)*D+I/p.Rb_ohm(j)]\[w.u;zeros(63,1)];
        phi(:,j)=h(1:63);
    end
    for target=1:4
        v=w.referencias{target}; x=v.x; on=bits(target,:);
        icp=p.Cp_F*D*x(:,2);
        residual=[w.e-v.u-cfg.Fuente.R_ohm*v.is-cfg.Fuente.L_H*D*v.is; ...
            v.u-x(:,2)-p.Rp_ohm*icp;v.is-il-icp-x(:,3)-x(:,5)];
        if cfg.Fuente.L_H>0, residual=[residual;x(:,1)-v.is]; else, residual=[residual;x(:,1)]; end %#ok<AGROW>
        for j=1:2
            ii=x(:,2*j+1); vc=x(:,2*j+2);
            if on(j)
                residual=[residual;v.u-p.R_ohm(j)*ii-vc-p.L_H(j)*D*ii; ...
                    ii-p.C_F(j)*D*vc-vc/p.Rb_ohm(j)]; %#ok<AGROW>
            else, residual=[residual;ii;vc]; %#ok<AGROW>
            end
        end
        q.Residuo_periodico=max(q.Residuo_periodico,max(abs(residual)));
        costs(target)=mean(v.is.^2); frozen(target)=mean((il+icpFixed+phi*on').^2);
    end
    q.Error_costes_A2=max([q.Error_costes_A2,abs(costs-d.CostesRed(k,:)),abs(frozen-d.CostesBornes(k,:))]);
    actual=d.EstadoActual(k)+1; chosen=d.EstadoSeleccionado(k)+1;
    control=costs; if c.ModeloCoste=="bornes", control=frozen; end
    q.Error_costes_A2=max(q.Error_costes_A2,max(abs(control-d{k,{'Coste00','Coste01','Coste10','Coste11'}})));
    reachable=sum(abs(bits-bits(actual,:)),2)<=1;
    assert(isequal(reachable,logical(d.Alcanzables(k,:)')) && d.Factibles(k,actual), ...
        'TSCRed:Vecindad','Vecindad o permanencia en el estado actual incorrecta.');
    eligible=reachable & logical(d.Factibles(k,:)');
    work=control; work(~eligible)=Inf;
    assert(eligible(chosen) && abs(control(chosen)-min(work))<1e-8, ...
        'TSCRed:Seleccion','No se selecciona un minimo entre vecinos admisibles.');
    candidate=d.EstadoCandidato(k)+1;
    assert(abs(control(candidate)-min(control))<1e-8,'TSCRed:Candidato','Candidato global incorrecto.');
    gain=control(actual)-control(chosen);
    assert(abs(gain-d.Ganancia_A2(k))<1e-8,'TSCRed:Ganancia','Ganancia incorrecta.');
    changed=~d.Pendiente(k) && d.Permanencia_s(k)>=p.permanenciaMin_s-1e-12 && ...
        gain>p.histeresis_A2 && actual~=chosen;
    assert(d.Cambio(k)==changed,'TSCRed:Retencion','No respeta histeresis/permanencia/pendiente.');
    [distance,it]=min(abs(a.t-td)); assert(distance<1e-9);
    xs=[a.iFuente(it),a.vCp(it),a.iRamas(it,1),a.vCs(it,1),a.iRamas(it,2),a.vCs(it,2)];
    if cfg.Fuente.L_H==0, xs(1)=0; end
    q.Error_estado_actual=max(q.Error_estado_actual,max(abs(xs-d.EstadoMedido(k,:))));
    if changed
        assert(d.EstadoSolicitado(k)==chosen-1 && nnz(bits(actual,:)~=bits(chosen,:))==1, ...
            'TSCRed:Orden','La orden no corresponde al estado seleccionado.');
        j=find(bits(chosen,:)>bits(actual,:));
        if ~isempty(j)
            tau=d.(sprintf('ProgramacionON%d_s',j))(k)-td;
            assert(tau>=-1e-10 && tau<=p.T+1e-10,'TSCRed:Plan','Instante fuera del horizonte.');
            old=w.referencias{actual}.x; new=w.referencias{chosen}.x;
            % Propagacion independiente en coordenadas fisicas (el control
            % usa matrices precalculadas en coordenadas de energia).
            xp=interpolar_periodico(old,p.omega*(td+tau))+ ...
                (expm(c.topologias{actual}.A*tau)*(d.EstadoMedido(k,:)-old(1,:))')';
            target=interpolar_periodico(new,p.omega*(td+tau));
            en=.5*sum(c.pesos.*(xp-target).^2);
            mv=referenceMaximum(new(:,[4 6]),cfg); mi=referenceMaximum(new(:,[3 5]),cfg);
            [vb,ib]=bounds(xp,en,mv,mi,bits(chosen,:),cfg);
            assert(all(vb<=(1-c.Margen)*c.MaxVC_V+1e-6) && ...
                all(ib<=(1-c.Margen)*c.MaxI_A+1e-7), ...
                'TSCRed:PlanLimites','La orden no respeta la cota prevista.');
            assert(abs(en-d.EnergiasCandidatas_J(k,chosen))<1e-6 && ...
                abs(max(vb)-d.CotasVC_V(k,chosen))<1e-5 && ...
                abs(max(ib)-d.CotasI_A(k,chosen))<1e-6, ...
                'TSCRed:PlanMagnitudes','Cotas/energia planificadas inconsistentes.');
        end
    end
end
assert(q.Error_captura<1e-7 && q.Error_estado_actual<1e-6,'TSCRed:Captura','Medidas falsas o desplazadas.');
assert(q.Residuo_periodico<1e-6,'TSCRed:Modelo','Las referencias no satisfacen el circuito acoplado.');
assert(q.Error_costes_A2<1e-7,'TSCRed:Costes','Los costes no son las normas RMS declaradas.');
contextRows=rows(r.contextos_admision,385); guards=struct([]);
for j=1:2
    g=rows(r.guardias{j},15); indices=find(diff([0;g(:,1)])>0);
    events=r.eventos(r.eventos.Rama==j & ismember(r.eventos.Tipo,["ON","rechazo_ON"]),:);
    assert(numel(indices)==height(events),'TSCRed:GuardiaRegistro','Faltan comprobaciones de ON.');
    for n=1:numel(indices)
        v=g(indices(n),:); ta=v(2); accepted=logical(v(3));
        assert(abs(ta-events.Tiempo_s(n))<1e-9 && accepted==(events.Tipo(n)=="ON"), ...
            'TSCRed:GuardiaEvento','La admision no corresponde al evento real.');
        planned=d.(sprintf('ProgramacionON%d_s',j)); [distance,k]=min(abs(planned-ta));
        assert(distance<1e-8 && d.Tiempo_s(k)<=ta+1e-10,'TSCRed:GuardiaCausal','ON sin plan antecedente.');
        [distance,ic]=min(abs(r.contextos_admision.Time-d.Tiempo_s(k))); assert(distance<1e-9);
        ctx=contextRows(ic,:)'; target=ctx(2)+1;
        assert(ctx(1)==j && ctx(3)<=ta+1e-10 && target==d.EstadoSeleccionado(k)+1, ...
            'TSCRed:Contexto','Contexto de revalidacion incorrecto.');
        reference=reshape(ctx(4:381),63,6); declared=r.ventanas_decision{k}.referencias{target}.x;
        assert(max(abs(reference-declared),[],'all')<1e-9,'TSCRed:Contexto','Referencia sustituida.');
        xp=interpolar_periodico(reference,p.omega*ta); x=v(10:15);
        [~,ix]=min(abs(a.t-ta)); actualX=[a.iFuente(ix),a.vCp(ix),a.iRamas(ix,1),a.vCs(ix,1),a.iRamas(ix,2),a.vCs(ix,2)];
        if cfg.Fuente.L_H==0, actualX(1)=0; end
        assert(max(abs(x-actualX))<1e-5,'TSCRed:EstadoEvento','Estado de revalidacion falso.');
        energy=.5*sum(c.pesos.*(x-xp).^2);
        mv=referenceMaximum(reference(:,[4 6]),cfg); mi=referenceMaximum(reference(:,[3 5]),cfg);
        assert(max(abs([mv mi]-ctx(382:385)'))<1e-8,'TSCRed:CotaReferencia','Maximos de referencia falsos.');
        [vb,ib]=bounds(x,energy,mv,mi,bits(target,:),cfg);
        q.Error_energia_evento_J=max(q.Error_energia_evento_J,abs(energy-v(5)));
        q.Error_cota_evento=max([q.Error_cota_evento,abs([vb ib]-v(6:9))]);
        expected=all(vb<=c.MaxVC_V) && all(ib<=c.MaxI_A) && ...
            abs(events.DiDt_A_s(n))<=p.limiteDiDt_A_s+1e-7*max(1,p.limiteDiDt_A_s);
        assert(accepted==expected,'TSCRed:Admision','Se admite/rechaza un ON contra las restricciones.');
        q.ON_admitidos=q.ON_admitidos+accepted; q.ON_rechazados=q.ON_rechazados+~accepted;
        guards=[guards;struct('Rama',j,'Tiempo_s',ta,'Admitido',accepted,'Razon',v(4), ...
            'Energia_J',energy,'CotaVC1_V',vb(1),'CotaVC2_V',vb(2),'CotaI1_A',ib(1),'CotaI2_A',ib(2))]; %#ok<AGROW>
    end
end
assert(q.Error_cota_evento<1e-6 && q.Error_energia_evento_J<1e-7,'TSCRed:GuardiaCota','Registro de cotas incoherente.');
q.PicosVC_V=max(abs(a.vCs),[],1); q.PicosI_A=max(abs(a.iRamas),[],1);
q.Limites_cumplidos_en_traza=all(q.PicosVC_V<=c.MaxVC_V) && all(q.PicosI_A<=c.MaxI_A);
q.ExcesoVC_V=max(0,q.PicosVC_V-c.MaxVC_V); q.ExcesoI_A=max(0,q.PicosI_A-c.MaxI_A);
if opciones.InspeccionarModelo, inspect(r); end
report=struct('passed',true,'fisica',base,'control',q,'admisiones',struct2table(guards), ...
    'alcance',"Verifica prediccion/admisiones y mide picos observados. La cota es condicional a excitaciones repetidas; no garantiza cambios futuros arbitrarios.");
end

function mv=referenceMaximum(v,cfg)
N=cfg.ControlRed.MuestrasCota; phi=(0:N-1)'*2*pi/N;
mv=max(abs(interpolar_periodico(v,phi)),[],1)+pi/N*sqrt(sum((derivador_periodico(63)*v).^2,1));
end
function [vb,ib]=bounds(x,en,mv,mi,bits,cfg)
vb=mv+sqrt(2*en./cfg.p.C_F); ib=mi+sqrt(2*en./cfg.p.L_H);
for j=1:2
    if ~bits(j), vb(j)=abs(x(2*j+2)); ib(j)=0; end
end
end
function v=rows(ts,width)
data=ts.Data; n=numel(ts.Time);
if size(data,1)==n && size(data,2)==width, v=reshape(data,n,width); else, v=reshape(data,width,n)'; end
end
function inspect(r)
load_system(r.modelo); [~,name]=fileparts(r.modelo); cleanup=onCleanup(@()close_system(name,0));
assert(strcmp(get_param([name,'/Supervisor_temporal'],'FunctionName'),'tsc_supervisor_red_sfun'));
sources=find_system(name,'BlockType','FromWorkspace');
assert(isscalar(sources) && strcmp(get_param(sources{1},'VariableName'),'tm_carga'),'TSCRed:Estructura','Cronograma externo detectado.');
assert(strcmp(get_param(name,'AlgebraicLoopMsg'),'error'),'TSCRed:Estructura','Bucles no auditados.');
ports=get_param([name,'/Medidas'],'PortConnectivity');
assert(strcmp(get_param(ports(1).SrcBlock,'Name'),'Tension_bornes') && ...
    strcmp(get_param(ports(8).SrcBlock,'Name'),'Tension_continua'),'TSCRed:Medidas','Faltan medidas de ambos lados de la red.');
for j=1:2
    ex=[name,sprintf('/Ejecutor_%d',j)];
    assert(strcmp(get_param([ex,'/Latch_evento/Actualizar'],'FunctionName'),'tsc_ejecutor_red_sfun'));
    assert(strcmp(get_param([ex,'/Cero_corriente'],'BlockType'),'HitCross'));
    assert(strcmp(get_param([name,sprintf('/LC_%d/Estado_VC',j)],'ExternalReset'),'none'));
end
end
