% Motor conservado de ejemplo_conmutacion_tsc.m; mismas ecuaciones y eventos.
function a=paper_maniobra(p,policy,label)
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
