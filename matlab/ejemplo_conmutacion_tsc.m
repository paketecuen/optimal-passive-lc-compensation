function result = ejemplo_conmutacion_tsc(opciones)
%EJEMPLO_CONMUTACION_TSC Maniobras de una LC fija con una valvula idealizada.
% Modelo hibrido temporal, no modelo de semiconductor ni ensayo Simulink.
% Se conservan los estados entre maniobras. La retirada de puerta espera
% el siguiente cero de corriente; el condensador conserva su tension.
% La otra rama Cp queda permanentemente conectada: no se conmuta C puro.
arguments
    opciones.Visible (1,1) string = "on"
    opciones.Exportar (1,1) logical = true
    opciones.PasoMaximo (1,1) double {mustBePositive} = 1/40000
end
c=caso_m3(63); f=familia_optima_m3(c.u,c.i,c.D,c.omega);
b=banco_familia_m3(f,27/13);
p=struct('L',b.Ls_H,'C',b.Cs_F,'Cp',b.Cp_F,'omega',c.omega, ...
    'T',2*pi/c.omega,'c',c,'request',[.021,.101,.161,.251], ...
    'end',.300,'h',opciones.PasoMaximo);
policies=["inmediato","valvula","estado","estado","valvula"];
labels=["Inmediato con perdidas","Minimo voltaje valvula", ...
    "Minima energia de desajuste","Ideal sin amortiguacion","Carga residual alta"];
runs=cell(5,1); summary=struct([]); events=table;
for k=1:5
    pk=p; pk.R=.05; pk.Rb=1000; pk.v0=0;
    if k==4, pk.R=0; pk.Rb=Inf; end
    if k==5, pk.v0=1.5*max(abs(tension(linspace(0,p.T,10001),p))); end
    runs{k}=simular(pk,policies(k),labels(k)); a=runs{k};
    summary=[summary;struct('Caso',labels(k),'R_ohm',pk.R,'Rb_ohm',pk.Rb, ...
        'VC_inicial_V',pk.v0,'I_serie_pico_A',max(abs(a.x(:,1))), ...
        'VC_pico_V',max(abs(a.x(:,2))), ...
        'V_valvula_bloqueada_max_V',max(abs(a.voltaje_valvula(~a.on))), ...
        'Error_balance_J',a.balance_error,'Error_apagado_A',a.off_error, ...
        'Salto_VC_V',a.voltage_jump)]; %#ok<AGROW>
    events=[events;a.events]; %#ok<AGROW>
end
summary=struct2table(summary);
% Pruebas fisicas: continuidad, corriente de extincion y balance energetico.
assert(all(summary.Error_balance_J<2e-8));
assert(all(summary.Error_apagado_A<1e-8));
assert(all(summary.Salto_VC_V<1e-12));
for k=1:5
    e=runs{k}.events;
    assert(all(e.Retardo_ms>=-1e-9));
    assert(all(e.Retardo_ms(e.Accion=="conectar")<=1000*p.T+1e-6));
end
assert(runs{5}.events.Valvula_min_posible_V(1)>0, ...
    'TSC:Residual','El caso cargado debe tener tension residual superior a la fuente.');
ideal=runs{4}; intervals=ideal.events.Tiempo_s;
constant_error=0;
for j=[1,3]
    idx=ideal.t>intervals(j)+p.h & ideal.t<intervals(j+1)-p.h;
    constant_error=max(constant_error,max(ideal.error_energy(idx))-min(ideal.error_energy(idx)));
end
assert(constant_error<2e-8,'TSC:Ideal','Un LC ideal debe conservar la energia libre.');
figures=gobjects(3,1);
figures(1)=figure('Visible',opciones.Visible,'Color','w','Position',[80,80,1250,850]);
tiledlayout(3,1,'TileSpacing','compact'); colors=lines(3);
nexttile; hold on
for k=1:3, plot(1000*runs{k}.t,runs{k}.x(:,1),'Color',colors(k,:)); end
ylabel('Corriente serie (A)'); legend(labels(1:3),'Location','best'); grid on
title('Dos conexiones y dos retiradas: la LC mantiene sus estados');
nexttile; hold on
plot(1000*runs{3}.t,runs{3}.u,'k:','DisplayName','Tension de bornes');
for k=1:3, plot(1000*runs{k}.t,runs{k}.x(:,2),'Color',colors(k,:),'DisplayName',labels(k)); end
ylabel('Tension (V)'); grid on
nexttile; hold on
for k=1:3, stairs(1000*runs{k}.t,double(runs{k}.on)+.08*(k-2),'Color',colors(k,:)); end
for t=p.request, xline(1000*t,'k:'); end
ylim([-.3,1.3]); ylabel('ON: trazas desplazadas'); xlabel('Tiempo (ms)'); grid on
title('Lineas verticales: orden externa; apagar exige esperar corriente cero');
figures(2)=figure('Visible',opciones.Visible,'Color','w','Position',[110,110,1250,850]);
tiledlayout(2,2,'TileSpacing','compact');
nexttile; hold on
for k=1:3
    a=runs{k}; plot(sqrt(p.L)*a.x(:,1),sqrt(p.C)*a.x(:,2),'Color',colors(k,:));
end
xlabel('sqrt(L) i'); ylabel('sqrt(C) v_C'); grid on; axis equal
title('Estado en el marco {P_+,N_-}');
nexttile; hold on
for k=[3,4]
    a=runs{k}; e=a.error_energy; e(~a.on)=NaN;
    plot(1000*a.t,1000*e,'DisplayName',labels(k));
end
xlabel('Tiempo (ms)'); ylabel('Energia libre respecto a forzada (mJ)');
legend('Location','best'); grid on
title('Las perdidas amortiguan; el ideal conserva la oscilacion libre');
nexttile; hold on
for k=1:3
    a=runs{k}; energy=.5*p.L*a.x(:,1).^2+.5*p.C*a.x(:,2).^2;
    balance=energy-energy(1)-a.x(:,3)+a.x(:,4);
    plot(1000*a.t,1e12*balance,'Color',colors(k,:));
end
xlabel('Tiempo (ms)'); ylabel('Residuo balance (pJ)'); grid on
title('E(t)-E(0) = integral(u i - perdidas) dt');
nexttile; a=runs{5}; hold on
plot(1000*a.t,a.u,'k:'); plot(1000*a.t,a.x(:,2));
xlim([0,55]); xlabel('Tiempo (ms)'); ylabel('Tension (V)'); grid on
legend('Bornes','Condensador residual','Location','best');
title('Si la carga residual excede la fuente no existe cruce de valvula');
figures(3)=figure('Visible',opciones.Visible,'Color','w','Position',[140,140,1250,550]);
tiledlayout(1,2,'TileSpacing','compact');
for k=1:2
    nexttile; a=runs{k+1}; scan=a.scans{1}; yyaxis left
    plot(1000*(scan.t-p.request(1)),1000*scan.energy); ylabel('Energia desajuste (mJ)');
    yyaxis right; plot(1000*(scan.t-p.request(1)),abs(scan.valve)); ylabel('|u-v_C| (V)');
    xline(1000*(a.events.Tiempo_s(1)-p.request(1)),'k--');
    xlabel('Espera desde orden (ms)'); grid on; title(labels(k+1));
end
result=struct('parametros',p,'banco',b,'resumen',summary,'eventos',events, ...
    'registros',{runs},'figuras',figures,'error_energia_ideal_J',constant_error, ...
    'alcance','Valvula idealizada; no semiconductor, Simulink ni hardware');
if opciones.Exportar
    folder=fullfile(fileparts(mfilename('fullpath')),'resultados','conmutacion_tsc');
    if ~isfolder(folder), mkdir(folder); end
    writetable(summary,fullfile(folder,'resumen.csv'));
    writetable(events,fullfile(folder,'eventos.csv'));
    names=["01_maniobras","02_estados_energia","03_criterios_disparo"];
    for k=1:3
        exportgraphics(figures(k),fullfile(folder,names(k)+".png"),'Resolution',150);
        savefig(figures(k),fullfile(folder,names(k)+".fig"));
    end
    saved=rmfield(result,'figuras'); save(fullfile(folder,'conmutacion_tsc.mat'),'saved');
    report=struct('passed',true,'modelo',result.alcance, ...
        'balance_max_J',max(summary.Error_balance_J), ...
        'apagado_max_A',max(summary.Error_apagado_A), ...
        'salto_condensador_max_V',max(summary.Salto_VC_V), ...
        'energia_ideal_variacion_max_J',constant_error);
    fid=fopen(fullfile(folder,'verificacion.json'),'w');
    cleanup=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
disp(summary); disp(events);
end

function a=simular(p,policy,label)
N=numel(p.c.u); Dt=p.omega*p.c.D; I=eye(N);
reference=[p.L*Dt+p.R*I,I;-I,p.C*Dt+I/p.Rb]\[p.c.u;zeros(N,1)];
ref=[reference(1:N),reference(N+1:end)];
t=0; x=[0;p.v0;0;0]; tt=0; xx=x'; modes=false; events=struct([]);
scans=cell(2,1); offerror=0; jump=0;
for n=1:4
    request=p.request(n); action=mod(n,2)==1;
    if action
        advance(request,false,false);
        [firing,scan]=schedule(t,x(2),p,ref,policy); scans{(n+1)/2}=scan;
        advance(firing,false,false);
        val=tension(t,p)-x(2); xp=interpolar_periodico(ref,p.omega*t);
        E=.5*p.L*(x(1)-xp(1))^2+.5*p.C*(x(2)-xp(2))^2;
        append_event("conectar",request,E,val,min(abs(scan.valve)));
    else
        advance(request,true,false);
        if abs(x(1))>1e-12, advance(request+2*p.T,true,true); end
        offerror=max(offerror,abs(x(1)));
        assert(abs(x(1))<1e-8,'TSC:Extincion','No se permite apertura a corriente no nula.');
        before=x(2); x(1)=0; jump=max(jump,abs(x(2)-before));
        append_event("bloquear",request,NaN,tension(t,p)-x(2),NaN);
    end
end
advance(p.end,false,false);
u=tension(tt,p); reference=interpolar_periodico(ref,p.omega*tt);
Ee=.5*p.L*(xx(:,1)-reference(:,1)).^2+.5*p.C*(xx(:,2)-reference(:,2)).^2;
energy=.5*p.L*xx(:,1).^2+.5*p.C*xx(:,2).^2;
vv=u-xx(:,2); vv(modes)=0;
a=struct('t',tt,'x',xx,'on',modes,'u',u,'voltaje_valvula',vv, ...
    'reference',reference,'error_energy',Ee,'scans',{scans}, ...
    'events',struct2table(events),'off_error',offerror,'voltage_jump',jump, ...
    'balance_error',max(abs(energy-energy(1)-xx(:,3)+xx(:,4))));
a.corriente_fuente=interpolar_periodico(p.c.i,p.omega*tt)+p.Cp*derivada_tension(tt,p)+xx(:,1);

    function advance(target,on,stop)
        if target<=t+1e-14, return; end
        options=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',p.h);
        if stop, options=odeset(options,'Events',@zero_current); end
        sol=ode45(@(ta,xa)dynamics(ta,xa,p,on),[t,target],x,options);
        tend=sol.x(end);
        if stop, assert(isfield(sol,'xe') && ~isempty(sol.xe),'TSC:Cero','No se encontro cero de corriente.'); end
        times=(t+p.h:p.h:tend)'; times=[times;tend];
        times=unique(times); values=deval(sol,times)';
        tt=[tt;times]; xx=[xx;values]; modes=[modes;repmat(on,numel(times),1)];
        t=tend; x=sol.y(:,end);
    end
    function append_event(action,request,E,val,vmin)
        didt=NaN; if action=="conectar", didt=val/p.L; end
        events=[events;struct('Caso',label,'Accion',action,'Orden_s',request, ...
            'Tiempo_s',t,'Retardo_ms',1000*(t-request),'I_A',x(1),'VC_V',x(2), ...
            'Valvula_V',val,'dI_dt_inicial_A_s',didt,'dU_dt_V_s',derivada_tension(t,p), ...
            'Desajuste_mJ',1000*E,'Valvula_min_posible_V',vmin)];
    end
end

function dx=dynamics(t,x,p,on)
if on
    dx=[(tension(t,p)-p.R*x(1)-x(2))/p.L; ...
        (x(1)-x(2)/p.Rb)/p.C;tension(t,p)*x(1);p.R*x(1)^2+x(2)^2/p.Rb];
else
    dx=[0;-x(2)/(p.Rb*p.C);0;x(2)^2/p.Rb];
end
end

function [value,terminal,direction]=zero_current(~,x)
value=x(1); terminal=1; direction=0;
end

function [firing,scan]=schedule(t0,v0,p,reference,policy)
tt=linspace(t0,t0+p.T,2001)'; vc=v0*exp(-(tt-t0)/(p.Rb*p.C));
ref=interpolar_periodico(reference,p.omega*tt);
energy=.5*p.L*ref(:,1).^2+.5*p.C*(vc-ref(:,2)).^2;
valve=tension(tt,p)-vc;
scan=struct('t',tt,'energy',energy,'valve',valve);
if policy=="inmediato", firing=t0; return; end
if policy=="estado", cost=energy; else, cost=valve.^2; end
candidates=find(cost(2:end-1)<=cost(1:end-2) & cost(2:end-1)<=cost(3:end))+1;
ts=[tt(1);tt(end)]; cs=[cost(1);cost(end)];
for j=candidates'
    [tj,cj]=fminbnd(@objective,tt(j-1),tt(j+1),optimset('TolX',1e-13));
    ts(end+1,1)=tj; cs(end+1,1)=cj; %#ok<AGROW>
end
[ts,order]=sort(ts); cs=cs(order);
best=min(cs); idx=find(cs<=best+1e-13*max(1,max(cs)),1);
firing=ts(idx);
    function y=objective(t)
        v=v0*exp(-(t-t0)/(p.Rb*p.C));
        if policy=="valvula", y=(tension(t,p)-v)^2; return; end
        r=interpolar_periodico(reference,p.omega*t);
        y=.5*p.L*r(1)^2+.5*p.C*(v-r(2))^2;
    end
end

function u=tension(t,p)
a=cos(p.omega*t); u=sqrt(2)*(4*a.^3+2*a.^2-2*a-1);
end

function du=derivada_tension(t,p)
a=cos(p.omega*t); s=sin(p.omega*t);
du=p.omega*sqrt(2)*(-12*a.^2.*s-4*a.*s+2*s);
end
