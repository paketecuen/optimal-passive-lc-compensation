function [b,status,extra] = paper_disenar(c,method)
%PAPER_DISENAR Reimplementaciones declaradas sobre datos comunes del paper.
N=c.N; D=c.D; v=D*c.u; j=c.J*c.u; extra=struct;
Dr=c.Dr; Ir=eye(size(Dr)); ur=c.V'*c.u/N;
b=struct('ramas',zeros(0,3)); status="realizable";
switch method
    case "Sin compensar"
        return
    case "C Kusters-Moore"
        a=max(0,-dot(c.i,v)/dot(v,v));
        if a>0, b.ramas=[0 a/c.omega 0]; end
        if c.Rg>0 || c.Lg>0, status="formula de tension fijada aplicada a u inicial; no es optimo con red"; end
    case "C optimizada circuito"
        if c.Rg==0 && c.Lg==0
            [b,status,extra]=paper_disenar(c,"C Kusters-Moore"); return
        end
        adm=norm(c.i)/norm(c.u); grid=linspace(log(1e-6*adm),log(1e3*adm),180);
        vals=arrayfun(@(x)fixedcost(x,0),grid); candidates=find(vals(2:end-1)<=vals(1:end-2) & vals(2:end-1)<=vals(3:end))+1;
        best=norm(c.i)^2/N; C0=0;
        for a=candidates
            [x,value]=fminbnd(@(x)fixedcost(x,0),grid(a-1),grid(a+1),optimset('TolX',1e-11));
            if value<best, best=value; C0=exp(x); end
        end
        if C0>0, b.ramas=[0 C0/c.omega 0]; end
        status="C optimizada numericamente con circuito acoplado; intervalo finito declarado";
    case "Paralelo Page"
        fit=cono_temporal([v j],c.i); a=fit.alpha;
        if a(1)>0, b.ramas(end+1,:)=[0 a(1)/c.omega 0]; end
        if a(2)>0, b.ramas(end+1,:)=[1/(c.omega*a(2)) 0 0]; end
        if c.Rg>0 || c.Lg>0, status="cono de tension fijada aplicado a u inicial; no es optimo con red"; end
    case "Formula m2"
        if numel(c.n)~=2, status="no aplica: dimension distinta de dos"; return; end
        sol=resolver_lc_m2(c.u,c.i,D,c.omega); extra=rmfield(sol,{'D','u','i'});
        if ~sol.serie_admisible, status="q* no realizable por una LC positiva finita"; return; end
        b.ramas=[sol.L_H sol.C_F 0];
    case "Busqueda directa LC"
        adm=norm(c.i)/norm(c.u); box=log(1e4); best=Inf;
        poles=c.n.^2; targets=[poles(1)/2,(poles(1:end-1)+poles(2:end))/2,2*poles(end)];
        options=optimset('Display','off','MaxIter',700,'MaxFunEvals',1500,'TolX',1e-9,'TolFun',1e-11);
        for L0=[.03 .3 3 30 300]
            for z=targets
                x0=log([L0,1/(z*L0)]); x0=max(-.95*box,min(.95*box,x0));
                [x,value]=fminsearch(@cost,x0,options);
                if value<best, best=value; chosen=x; end
            end
        end
        b.ramas=[exp(chosen(1))/(adm*c.omega),exp(chosen(2))*adm/c.omega,0];
        status="busqueda local multiinicio en caja L/C de ocho decadas; no certificado global";
        extra=struct('log_parametros',chosen,'log_limite',box,'toca_caja',any(abs(chosen)>.99*box));
        if ismember(c.id,["Jeltsema2014","Czarnecki1987_ideal"])
            status="aproximacion finita al infimo en la frontera; esfuerzo no acotado en el limite";
        elseif c.id=="Czarnecki1987_red"
            status="aproximacion a frontera resonante; esfuerzo no acotado; globalidad con red no certificada";
        elseif startsWith(c.id,"Czarnecki1981")
            status="minimo global LC positivo finito; auditoria exacta independiente de candidatos y fronteras";
        end
    case "Busqueda LC restringida"
        [b,status,extra]=paper_lc_restringida(c);
    case "Czarnecki L fija"
        adm=norm(c.i)/norm(c.u); Lgrid=[.01 .03 .1 .3 1 2 3 10]/adm;
        if startsWith(c.id,"Czarnecki1987"), Lgrid=[.06 .10 .20]*c.omega; end
        % Resuelve el problema de C con L fijada por intervalos entre polos.
        % Para fuente ideal es la estacionariedad de ec.8 de 1987; con red
        % se minimiza directamente el circuito acoplado, no se usa ec.9.
        best=Inf; gridResults=zeros(numel(Lgrid),3);
        for ell=1:numel(Lgrid)
            L0=Lgrid(ell); limits=sort([1e-6*adm,1./(L0*c.n.^2),1e4*adm]); localBest=Inf;
            for k=1:numel(limits)-1
                lo=log(limits(k)*(1+1e-7)); hi=log(limits(k+1)*(1-1e-7));
                xx=linspace(lo,hi,65); yy=arrayfun(@(x)fixedcost(x,L0),xx);
                [value,idx]=min(yy); xc=xx(idx);
                centers=find(yy(2:end-1)<=yy(1:end-2) & yy(2:end-1)<=yy(3:end))+1;
                for center=centers
                    [x,vv]=fminbnd(@(x)fixedcost(x,L0),xx(center-1),xx(center+1),optimset('TolX',1e-11));
                    if vv<value, value=vv; xc=x; end
                end
                if value<localBest, localBest=value; localC=exp(xc); end
            end
            gridResults(ell,:)=[L0/c.omega,localC/c.omega,sqrt(localBest)];
            if localBest<best, best=localBest; C0=localC; chosenL=L0; end
        end
        b.ramas=[chosenL/c.omega,C0/c.omega,0];
        extra=struct('L_fijada_H',chosenL/c.omega,'barrido_L_C_I',gridResults);
        status="mejor de varios L fijados; barrido escalar por intervalos; no optimo de L libre ni replica de iteracion (9)";
    case {"Bezout Cauer","Bezout Foster"}
        extra=paper_sintesis_bezout(c);
        if method=="Bezout Cauer", b=extra.cauer; else, b=extra.foster; end
    case {"Foster polos libres","Jeon polos previos"}
        % B_n se obtiene de la proyeccion temporal. La politica de polos
        % utiliza las frecuencias declaradas, como el metodo publicado.
        B=zeros(size(c.n));
        for k=1:numel(c.n)
            co=cos(c.n(k)*c.tau); si=sin(c.n(k)*c.tau);
            uc=2*dot(c.u,co)/N; us=-2*dot(c.u,si)/N;
            qc=2*dot(c.qstar,co)/N; qs=-2*dot(c.qstar,si)/N;
            B(k)=imag((qc+1i*qs)/(uc+1i*us));
        end
        offset=.2; if method=="Jeon polos previos", offset=.05; end
        for attempt=1:15
            mu=(c.n.*(1+offset*sign(B))).^2;
            if method=="Foster polos libres" && startsWith(c.id,"Czarnecki1987"), mu=[3.3 5.9 10.9].^2; end
            Phi=zeros(N,numel(mu));
            for k=1:numel(mu), Phi(:,k)=c.V*((Dr*Dr+mu(k)*Ir)\(Dr*ur)); end
            alpha=Phi\c.qstar;
            if all(alpha>0), break; end
            offset=offset/2;
        end
        if any(alpha<=0), status="polos elegidos no dan residuos positivos"; return; end
        b.ramas=[1./(c.omega*alpha),alpha./(c.omega*mu'),zeros(numel(mu),1)];
        extra=struct('mu',mu,'offset',offset,'intentos',attempt);
    case "LP polos fijos"
        % Familia con nucleo positivo y simplex, como SEGAN; se declara la
        % eleccion de polos y de cara (minimizar C paralelo) para reproducirla.
        mu=((c.n(1:end-1)+c.n(2:end))/2).^2;
        Phi=[v j];
        for k=1:numel(mu), Phi=[Phi,c.V*((Dr*Dr+mu(k)*Ir)\(Dr*ur))]; end %#ok<AGROW>
        [Q,~]=qr(Phi(:,1:numel(c.n)),0);
        scale=vecnorm(Phi); Aeq=Q'*(Phi./scale); rhs=Q'*c.qstar;
        objective=zeros(size(Phi,2),1); objective(1)=1/scale(1);
        [scaled,~,flag]=linprog(objective,[],[],Aeq,rhs,zeros(size(objective)),[], ...
            optimoptions('linprog','Display','none'));
        assert(flag>0,'Paper:LP','No se encontro solucion positiva.'); alpha=scaled./scale';
        if alpha(1)>1e-10, b.ramas(end+1,:)=[0 alpha(1)/c.omega 0]; end
        if alpha(2)>1e-10, b.ramas(end+1,:)=[1/(c.omega*alpha(2)) 0 0]; end
        for k=1:numel(mu)
            if alpha(k+2)>1e-10, b.ramas(end+1,:)=[1/(c.omega*alpha(k+2)) alpha(k+2)/(c.omega*mu(k)) 0]; end
        end
        extra=struct('mu',mu,'coeficientes',alpha,'residuo_A',norm(Phi*alpha-c.qstar)/sqrt(N));
    otherwise
        error('Paper:Metodo','Metodo no reconocido.');
end
    function value=cost(x)
        violation=max(0,abs(x)-box);
        if any(violation>0), value=1e3+sum(violation.^2); return; end
        bank=struct('ramas',[exp(x(1))/(adm*c.omega),exp(x(2))*adm/c.omega,0]);
        rr=paper_evaluar_red(c,bank,Trazas=false); value=rr.Is_RMS_A^2/(norm(c.i)^2/N);
        if ~isfinite(value), value=1e10; end
    end
    function value=fixedcost(x,Lbar)
        rr=paper_evaluar_red(c,struct('ramas',[Lbar/c.omega exp(x)/c.omega 0]),Trazas=false);
        value=rr.Is_RMS_A^2;
    end
end
