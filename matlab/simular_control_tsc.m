function r = simular_control_tsc(opciones)
%SIMULAR_CONTROL_TSC Supervisor causal y ejecutor hibrido de dos LC fijas.
% Planta temporal continua con fuente periodica impuesta. La carga a(t)*i0(t)
% es una corriente sintetica inyectada: no modela una carga pasiva variable.
% Valvula bidireccional ideal regateada, no semiconductor comercial.
% Cp tiene ESR y permanece conectado, inicialmente en respuesta periodica.
% Todas las estrategias comienzan con ambas LC bloqueadas y conservan VC.
arguments
    opciones.Modo (1,1) string {mustBeMember(opciones.Modo,["ga","convencional","fijo"])} = "ga"
    opciones.PNominal_W (1,1) double {mustBePositive,mustBeFinite} = 1000
    opciones.Duracion_s (1,1) double {mustBePositive,mustBeFinite} = .8
    opciones.PasoSalida_s (1,1) double {mustBePositive,mustBeFinite} = 1e-4
    opciones.Permanencia_s (1,1) double {mustBeNonnegative,mustBeFinite} = .06
    opciones.Histeresis_A2 (1,1) double {mustBeNonnegative,mustBeFinite} = .025
    opciones.MaxDiDt_A_s (1,1) double {mustBeNonnegative} = 2500
    opciones.MargenDiDt (1,1) double {mustBeNonnegative,mustBeFinite,mustBeLessThan(opciones.MargenDiDt,1)} = .01
    opciones.BitsFijos (1,2) double = [0 1]
    opciones.VCInicial_V (1,2) double {mustBeFinite} = [0 0]
    opciones.CoefTension (1,3) double {mustBeFinite} = [1 1 1]
    opciones.EscalaCarga (1,:) double {mustBePositive,mustBeFinite} = [.45 1.35 .55 1.15]
    opciones.InstantesCarga_s (1,:) double {mustBeNonnegative,mustBeFinite} = [0 .2 .4 .6]
    opciones.Transicion_s (1,1) double {mustBePositive,mustBeFinite} = .06
end
assert(all(ismember(opciones.BitsFijos,[0 1])), 'ControlTSC:Bits','BitsFijos debe contener ceros o unos.');
assert(numel(opciones.EscalaCarga)==numel(opciones.InstantesCarga_s) && ...
    opciones.InstantesCarga_s(1)==0 && all(diff(opciones.InstantesCarga_s)>=opciones.Transicion_s), ...
    'ControlTSC:Carga','Instantes de carga ordenados, primer instante cero y rampas sin solapar.');
c=caso_m3(63); scaled=dimensionar_referencia_230V(opciones.PNominal_W);
assert(opciones.CoefTension(1)>0,'ControlTSC:Tension','La componente fundamental debe ser positiva.');
co=cos(c.tau); coeff=opciones.CoefTension;
baseVoltage=sqrt(2)*(coeff(1)*co+coeff(2)*(2*co.^2-1)+coeff(3)*(4*co.^3-3*co));
primitive=[c.D,ones(numel(co),1);ones(1,numel(co)),0]\[baseVoltage;0];
baseCurrent=c.G*baseVoltage+c.C_normalizada*c.D*baseVoltage+c.gamma_L*primitive(1:end-1)+ ...
    c.alpha_serie*((c.D*c.D+c.z_serie*eye(numel(co)))\(c.D*baseVoltage));
b=scaled.banco; weight=[.5 1]; T=2*pi/c.omega;
p=struct('L_H',b.Ls_H./weight,'C_F',b.Cs_F*weight, ...
    'R_ohm',c.omega*b.Ls_H./weight/100+1./(c.omega*b.Cs_F*weight*500), ...
    'Rb_ohm',1000*scaled.escala_impedancia./weight, ...
    'Cp_F',b.Cp_F,'Rp_ohm',1/(c.omega*b.Cp_F*500), ...
    'omega',c.omega,'T',T,'lambda',230/norm(coeff),'k',scaled.escala_impedancia, ...
    'coefTension',coeff, ...
    'limiteDiDt_A_s',opciones.MaxDiDt_A_s, ...
    'limiteDiDtPlan_A_s',(1-opciones.MargenDiDt)*opciones.MaxDiDt_A_s, ...
    'margenDiDt',opciones.MargenDiDt,'permanenciaMin_s',opciones.Permanencia_s, ...
    'histeresis_A2',opciones.Histeresis_A2,'pasoMuestreo_s',opciones.PasoSalida_s, ...
    'ruidoSigma_A',0,'PNominal_W',opciones.PNominal_W,'opciones',opciones);
p.uReferencia=baseVoltage*p.lambda; p.iCargaReferencia=baseCurrent*p.lambda/p.k;
ref=referencias(p.uReferencia,c.D,p); p.cpReferencia=[ref.iCp ref.vCp];
p.x0=[0;opciones.VCInicial_V(1);0;opciones.VCInicial_V(2);0;0];
bits=[0 0;0 1;1 0;1 1]; t=0; x=p.x0; applied=false(1,2);
p.bits=bits;
requested=applied; pending=zeros(1,2); scheduled=[Inf Inf]; offDeadline=scheduled;
requestTime=[NaN NaN]; plannedSlope=[NaN NaN]; plannedEnergy=[NaN NaN];
lastConfirmed=-Inf; nextDecision=T;
tt=0; xx=x'; modes=applied; wants=requested; waits=pending;
events=struct([]); decisions=struct([]); windows=cell(0,1);
maxStep=min(opciones.PasoSalida_s,5e-5);
while t<opciones.Duracion_s-1e-13
    target=min([nextDecision,scheduled,offDeadline,opciones.Duracion_s]);
    hit=[];
    if target>t+1e-13
        opts=odeset('RelTol',2e-10,'AbsTol',1e-11,'MaxStep',maxStep);
        if any(pending==-1), opts=odeset(opts,'Events',@extincion); end
        sol=ode45(@dinamica,[t target],x,opts);
        tend=sol.x(end);
        if isfield(sol,'ie'), hit=sol.ie; end
        samples=((floor(t/opciones.PasoSalida_s+1e-8)+1)*opciones.PasoSalida_s:opciones.PasoSalida_s:tend)';
        samples=samples(samples>t+1e-13 & samples<tend-1e-13);
        samples=[samples;tend]; %#ok<AGROW>
        vals=deval(sol,samples)';
        tt=[tt;samples]; %#ok<AGROW>
        xx=[xx;vals]; %#ok<AGROW>
        modes=[modes;repmat(applied,numel(samples),1)]; %#ok<AGROW>
        wants=[wants;repmat(requested,numel(samples),1)]; %#ok<AGROW>
        waits=[waits;repmat(pending,numel(samples),1)]; %#ok<AGROW>
        t=tend; x=sol.y(:,end);
    end
    for j=1:2
        if pending(j)==-1 && (ismember(j,hit) || abs(x(2*j-1))<1e-9)
            before=x; x(2*j-1)=0; applied(j)=false; pending(j)=0; offDeadline(j)=Inf;
            lastConfirmed=t; appendEvent("OFF",j,before,x,NaN,NaN,true);
        end
    end
    for j=1:2
        if pending(j)==1 && scheduled(j)<=t+1e-12
            slope=(tension(t,p)-x(2*j))/p.L_H(j);
            predicted=interpolar_periodico(windows{end}.referenciaEstado{j},p.omega*t);
            free=.5*p.L_H(j)*predicted(1)^2+.5*p.C_F(j)*(x(2*j)-predicted(2))^2;
            if abs(slope)<=p.limiteDiDt_A_s+1e-7*max(1,p.limiteDiDt_A_s)
                applied(j)=true; pending(j)=0; lastConfirmed=t;
                appendEvent("ON",j,x,x,slope,free,true);
            else
                pending(j)=0; requested(j)=false;
                appendEvent("rechazo_ON",j,x,x,slope,free,false);
            end
            scheduled(j)=Inf;
        elseif pending(j)==-1 && offDeadline(j)<=t+1e-12
            % An extinction timeout does not authorize opening nonzero current.
            offDeadline(j)=Inf; appendEvent("timeout_OFF",j,x,x,NaN,NaN,false);
        end
    end
    if nextDecision<=t+1e-12 && t<opciones.Duracion_s-1e-13
        decidir(); nextDecision=nextDecision+T;
    end
    xx(end,:)=x'; modes(end,:)=applied; wants(end,:)=requested; waits(end,:)=pending;
end
u=tension(tt,p); icpvc=interpolar_periodico(p.cpReferencia,p.omega*tt);
ir=xx(:,[1 3]); vc=xx(:,[2 4]); loadScale=envolvente(tt,opciones);
il=loadScale.*interpolar_periodico(p.iCargaReferencia,p.omega*tt);
is=il+icpvc(:,1)+sum(ir,2);
loss=p.Rp_ohm*icpvc(:,1).^2+sum(ir.^2.*p.R_ohm+vc.^2./p.Rb_ohm,2);
energy=.5*p.Cp_F*icpvc(:,2).^2+.5*sum(ir.^2.*p.L_H+vc.^2.*p.C_F,2);
balance=energy-energy(1)-xx(:,5)+xx(:,6);
tr=struct('t',tt,'u',u,'iCarga',il,'iFuente',is,'iCp',icpvc(:,1),'iRamas',ir, ...
    'vCp',icpvc(:,2),'vCs',vc,'estadoSolicitado',wants,'estadoAplicado',modes, ...
    'pendiente',waits,'energia_J',energy,'integralEntrada_J',xx(:,5), ...
    'integralPerdidas_J',xx(:,6),'perdidas_W',loss,'escalaCarga',loadScale, ...
    'iCargaReconstruida',is-icpvc(:,1)-sum(ir,2));
ev=tablaEventos(events); de=tablaDecisiones(decisions);
summary=struct('Modo',opciones.Modo,'Duracion_s',tt(end),'U_RMS_V',sqrt(trapz(tt,u.^2)/tt(end)), ...
    'I_fuente_RMS_A',sqrt(trapz(tt,is.^2)/tt(end)), ...
    'I_sin_banco_RMS_A',sqrt(trapz(tt,il.^2)/tt(end)), ...
    'I_rama_pico_A',max(abs(ir),[],'all'),'VC_pico_V',max(abs(vc),[],'all'), ...
    'P_perdidas_media_W',xx(end,6)/tt(end),'Balance_max_J',max(abs(balance)), ...
    'Balance_trapz_final_J',energy(end)-energy(1)-trapz(tt,u.*(icpvc(:,1)+sum(ir,2))-loss), ...
    'Maniobras',sum(ev.Tipo=="ON" | ev.Tipo=="OFF"),'Rechazos',sum(ev.Tipo=="rechazo_ON"), ...
    'Pendientes_finales',nnz(pending),'Error_reconstruccion_A',max(abs(tr.iCargaReconstruida-il)), ...
    'Ramas_conmutables',2,'Ramas_fijas',1,'Componentes_LC',5);
r=struct('trazas',tr,'eventos',ev,'decisiones',de,'ventanas_decision',{windows}, ...
    'parametros',p,'referencias',ref,'resumen',summary, ...
    'alcance',"230 V y potencia nominal son semejanza de ensayo; carga inyectada, fuente fija, valvula ideal y Cp periodico. No certifica hardware ni una carga pasiva en transitorio.");

    function dx=dinamica(ta,xa)
        ua=tension(ta,p); cp=interpolar_periodico(p.cpReferencia,p.omega*ta);
        dx=zeros(6,1); diss=p.Rp_ohm*cp(1)^2; bankCurrent=cp(1);
        for branch=1:2
            ii=2*branch-1; vi=2*branch;
            if applied(branch)
                dx(ii)=(ua-p.R_ohm(branch)*xa(ii)-xa(vi))/p.L_H(branch);
                dx(vi)=(xa(ii)-xa(vi)/p.Rb_ohm(branch))/p.C_F(branch);
                bankCurrent=bankCurrent+xa(ii);
            else
                dx(vi)=-xa(vi)/(p.Rb_ohm(branch)*p.C_F(branch));
            end
            diss=diss+p.R_ohm(branch)*xa(ii)^2+xa(vi)^2/p.Rb_ohm(branch);
        end
        dx(5)=ua*bankCurrent; dx(6)=diss;
    end

    function [value,terminal,direction]=extincion(~,xa)
        value=ones(2,1); value(pending==-1)=xa(2*find(pending==-1)-1);
        terminal=ones(2,1); direction=zeros(2,1);
    end

    function appendEvent(kind,j,before,after,slope,free,admissible)
        events=[events;struct('Tipo',kind,'Rama',j,'Tiempo_s',t,'Orden_s',requestTime(j), ...
            'Retardo_s',t-requestTime(j),'I_pre_A',before(2*j-1),'I_post_A',after(2*j-1), ...
            'VC_pre_V',before(2*j),'VC_post_V',after(2*j), ...
            'DiDt_A_s',slope,'Desajuste_J',free,'Admisible',admissible, ...
            'DiDt_predicha_A_s',plannedSlope(j),'Desajuste_predicho_J',plannedEnergy(j))];
    end

    function decidir()
        N=numel(c.u); tw=t-T+(0:N-1)'*T/N;
        % Only observed samples at or before this decision are interpolated.
        measuredX=interp1(tt,xx(:,1:4),tw,'linear');
        measuredU=interp1(tt,tension(tt,p),tw,'linear');
        measuredCp=interp1(tt,interpolar_periodico(p.cpReferencia,p.omega*tt),tw,'linear');
        measuredLoad=interp1(tt,envolvente(tt,opciones).* ...
            interpolar_periodico(p.iCargaReferencia,p.omega*tt),tw,'linear');
        measuredBranches=measuredX(:,[1 3]);
        measuredSource=measuredLoad+measuredCp(:,1)+sum(measuredBranches,2);
        reconstructed=measuredSource-measuredCp(:,1)-sum(measuredBranches,2);
        forecast=referencias(measuredU,c.D,p); Phi=forecast.iRamas;
        v0=reconstructed+forecast.iCp; G=Phi'*Phi/N; h=Phi'*v0/N; constant=mean(v0.^2);
        costs=constant+2*(bits*h)'+sum((bits*G).*bits,2)';
        direct=mean((v0+Phi*bits').^2,1);
        cs=cos(p.omega*tw); sn=sin(p.omega*tw);
        ac=2*mean(measuredU.*cs); as=2*mean(measuredU.*sn);
        quadrature=sqrt(2)*(ac*sn-as*cs)/hypot(ac,as);
        if opciones.Modo=="convencional"
            controlCosts=mean((v0+Phi*bits').*quadrature,1).^2;
        else
            controlCosts=costs;
        end
        current=find(all(bits==applied,2)); best=min(controlCosts);
        ties=find(abs(controlCosts-best)<1e-10*max(1,best));
        [~,jTie]=min(sum(abs(bits(ties,:)-applied),2)); candidate=ties(jTie);
        if opciones.Modo=="fijo", candidate=find(all(bits==opciones.BitsFijos,2)); end
        age=t-lastConfirmed; wasPending=any(pending~=0);
        gain=controlCosts(current)-controlCosts(candidate);
        allow=~wasPending && age>=p.permanenciaMin_s-1e-12 && ...
            (gain>p.histeresis_A2 || opciones.Modo=="fijo");
        change=allow && candidate~=current;
        windows{end+1,1}=struct('t',tw,'u',measuredU,'iFuente',measuredSource, ...
            'iCp',measuredCp(:,1),'iRamas',measuredBranches,'iCargaReconstruida',reconstructed, ...
            'iCpReferencia',forecast.iCp,'phiRamas',Phi,'cuadraturaFundamental',quadrature, ...
            'G',G,'h',h,'c',constant,'costesRMS',costs,'costesControl',controlCosts, ...
            'referenciaEstado',{forecast.estado});
        if change
            requested=logical(bits(candidate,:));
            for branch=1:2
                if requested(branch)==applied(branch), continue; end
                requestTime(branch)=t;
                plannedSlope(branch)=NaN; plannedEnergy(branch)=NaN;
                if requested(branch)
                    appendEvent("solicitud_ON",branch,x,x,NaN,NaN,true);
                    [firing,slope,free]=programar(t,x(2*branch),branch,p,measuredU,forecast.estado{branch});
                    plannedSlope(branch)=slope; plannedEnergy(branch)=free;
                    if isfinite(firing)
                        pending(branch)=1; scheduled(branch)=firing;
                    else
                        requested(branch)=false;
                        currentSlope=(tension(t,p)-x(2*branch))/p.L_H(branch);
                        predicted=interpolar_periodico(forecast.estado{branch},p.omega*t);
                        currentFree=.5*p.L_H(branch)*predicted(1)^2+ ...
                            .5*p.C_F(branch)*(x(2*branch)-predicted(2))^2;
                        appendEvent("rechazo_ON",branch,x,x,currentSlope,currentFree,false);
                    end
                else
                    appendEvent("solicitud_OFF",branch,x,x,NaN,NaN,true);
                    pending(branch)=-1; offDeadline(branch)=t+2*T;
                end
            end
        end
        decisions=[decisions;struct('Tiempo_s',t,'InicioVentana_s',t-T,'FinVentana_s',t, ...
            'UltimoPuntoInterpolado_s',tw(end), ...
            'Coste00',controlCosts(1),'Coste01',controlCosts(2),'Coste10',controlCosts(3),'Coste11',controlCosts(4), ...
            'EstadoActual',current-1,'EstadoCandidato',candidate-1, ...
            'EstadoSolicitado',requested(1)*2+requested(2),'Pendiente',wasPending, ...
            'Permanencia_s',age,'Cambio',change,'Ganancia_A2',gain, ...
            'ErrorReconstruccion_A',max(abs(reconstructed-measuredLoad)), ...
            'ErrorGram_A2',max(abs(costs-direct)),'MaxMuestraUsada_s',tt(end))];
    end
end

function ref=referencias(u,D,p)
N=numel(u); Dt=p.omega*D; I=eye(N); states=cell(1,2); current=zeros(N,2);
for j=1:2
    solution=[p.L_H(j)*Dt+p.R_ohm(j)*I,I;-I,p.C_F(j)*Dt+I/p.Rb_ohm(j)]\[u;zeros(N,1)];
    states{j}=[solution(1:N),solution(N+1:end)]; current(:,j)=solution(1:N);
end
icp=(I+p.Rp_ohm*p.Cp_F*Dt)\(p.Cp_F*Dt*u);
ref=struct('iCp',icp,'vCp',u-p.Rp_ohm*icp,'iRamas',current,'estado',{states});
end

function [firing,slope,free]=programar(t0,vc0,j,p,uPast,reference)
% The last measured voltage period is repeated as a prediction. Load future
% is never queried. Each actual firing rechecks the instantaneous slope.
grid=linspace(t0,t0+p.T,1001)'; [energy,di]=objective(grid);
limitPlan=p.limiteDiDtPlan_A_s;
boundaries=grid([1 end]);
if isfinite(limitPlan)
    % Separate signed level crossings retain narrow admissible intervals and
    % zero-slope points, which a sign test on abs(di)-limit would lose.
    for level=unique([-limitPlan,limitPlan])
        signed=di-level; crossing=find(signed(1:end-1).*signed(2:end)<0);
        for k=crossing'
            boundaries(end+1,1)=fzero(@(ta)slopeOnly(ta)-level, ...
                [grid(k),grid(k+1)],optimset('TolX',1e-13)); %#ok<AGROW>
        end
        boundaries=[boundaries;grid(signed==0)]; %#ok<AGROW>
    end
    if limitPlan==0
        % A tangent zero need not change sign. This is a numerical search,
        % with the same tight feasibility tolerance as the other candidates.
        near=find(abs(di(2:end-1))<=abs(di(1:end-2)) & abs(di(2:end-1))<=abs(di(3:end)))+1;
        for k=near'
            boundaries(end+1,1)=fminbnd(@(ta)abs(slopeOnly(ta)), ...
                grid(k-1),grid(k+1),optimset('TolX',1e-13)); %#ok<AGROW>
        end
    end
end
candidate=boundaries; minima=find(energy(2:end-1)<=energy(1:end-2) & energy(2:end-1)<=energy(3:end))+1;
for k=minima'
    candidate(end+1,1)=fminbnd(@energyOnly,grid(k-1),grid(k+1),optimset('TolX',1e-13)); %#ok<AGROW>
end
[en,sl]=objective(candidate);
okay=abs(sl)<=limitPlan+1e-8*max(1,limitPlan);
candidate=candidate(okay); en=en(okay); sl=sl(okay);
if isempty(candidate)
    firing=Inf; [~,idx]=min(abs(di)); slope=di(idx); free=energy(idx); return
end
[free,idx]=min(en); firing=candidate(idx); slope=sl(idx);
    function [e,d]=objective(ta)
        vc=vc0*exp(-(ta-t0)/(p.Rb_ohm(j)*p.C_F(j)));
        % Decisions occur at integer periods; tau is the same phase as uPast.
        xp=interpolar_periodico(reference,p.omega*ta);
        up=interpolar_periodico(uPast,p.omega*ta);
        e=.5*p.L_H(j)*xp(:,1).^2+.5*p.C_F(j)*(vc-xp(:,2)).^2;
        d=(up-vc)/p.L_H(j);
    end
    function value=energyOnly(ta), value=objective(ta); end
    function value=slopeOnly(ta), [~,value]=objective(ta); end
end

function u=tension(t,p)
co=cos(p.omega*t); a=p.coefTension;
u=p.lambda*sqrt(2)*(a(1)*co+a(2)*(2*co.^2-1)+a(3)*(4*co.^3-3*co));
end

function a=envolvente(t,o)
a=repmat(o.EscalaCarga(1),size(t));
for j=2:numel(o.EscalaCarga)
    phase=max(0,min(1,(t-o.InstantesCarga_s(j))/o.Transicion_s));
    a=a+(o.EscalaCarga(j)-o.EscalaCarga(j-1))*(1-cos(pi*phase))/2;
end
end

function T=tablaEventos(s)
if isempty(s)
    T=table('Size',[0 14],'VariableTypes',{'string','double','double','double','double', ...
        'double','double','double','double','double','double','logical','double','double'}, ...
        'VariableNames',{'Tipo','Rama','Tiempo_s','Orden_s','Retardo_s','I_pre_A','I_post_A', ...
        'VC_pre_V','VC_post_V','DiDt_A_s','Desajuste_J','Admisible', ...
        'DiDt_predicha_A_s','Desajuste_predicho_J'});
else
    T=struct2table(s);
end
end

function T=tablaDecisiones(s)
if isempty(s)
    T=table;
else
    T=struct2table(s);
end
end
