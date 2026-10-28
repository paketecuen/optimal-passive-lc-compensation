function [bank,status,extra] = paper_lc_restringida(c)
%PAPER_LC_RESTRINGIDA LC positiva con cota de energia media del mismo caso.
% La cota es la menor energia de las CINCO realizaciones comparadas de orden
% libre; no un minimo universal de energia. Operadores temporales RMS.
references=["Bezout Cauer","Bezout Foster","Foster polos libres","Jeon polos previos","LP polos fijos"];
energies=zeros(size(references));
for k=1:numel(references)
    [b,~,~]=paper_disenar(c,references(k)); rr=paper_evaluar_red(c,b);
    energies(k)=rr.E_media_J;
end
[cap,idx]=min(energies); adm=norm(c.i)/norm(c.u); baseline=norm(c.i)^2/c.N;
D=c.Dr; I=eye(size(D)); Z=c.Rg*I+c.Lg*c.omega*D;
poles=c.n.^2; edges=[1e-4,poles,1e4*poles(end)];
best=1; chosen=[]; evaluations=0; flags=[]; starts=0;
options=optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',160,'MaxFunctionEvaluations',750,'OptimalityTolerance',1e-9, ...
    'StepTolerance',1e-11,'ConstraintTolerance',1e-9);
for region=1:numel(edges)-1
    lo=edges(region)*(1+1e-7); hi=edges(region+1)*(1-1e-7);
    lower=[log(lo),log(1e-7)]; upper=[log(hi),log(1e4)];
    if region==numel(edges)-1
        zs=exp(linspace(log(lo),log(hi),5)); zs=zs(2:4);
    else
        zs=lo+(.15:.35:.85)*(hi-lo);
    end
    for z=zs
        h=(D*D+z*I)\(D*(c.V'*c.u/c.N));
        alpha=max(adm*1e-4,-dot(c.V'*c.i/c.N,h)/dot(h,h));
        x0=[log(z),log(alpha/adm)]; x0=max(lower,min(upper,x0));
        [~,energy]=metrics(x0);
        while energy>cap*(1-1e-8) && x0(2)>lower(2)+log(2)
            x0(2)=x0(2)-log(2); [~,energy]=metrics(x0);
        end
        starts=starts+1;
        [x,value,flag]=fmincon(@objective,x0,[],[],[],[],lower,upper,@constraint,options);
        flags(end+1)=flag; %#ok<AGROW>
        % Un punto factible cuenta aun cuando la tolerancia de optimalidad
        % local no se alcanza. Se conserva tambien el inicio factible.
        consider(x,value); consider(x0,objective(x0));
    end
end
bank=struct('ramas',zeros(0,3));
if ~isempty(chosen)
    z=exp(chosen(1)); alpha=adm*exp(chosen(2));
    bank.ramas=[1/(c.omega*alpha),alpha/(c.omega*z),0];
end
status="LC con energia media limitada; busqueda local multiinicio por intervalos; sin certificado global";
extra=struct('E_limite_J',cap,'Referencia_energia',references(idx), ...
    'E_referencias_J',energies,'referencias',references,'z_limites',edges([1 end]), ...
    'alpha_sobre_adm_limites',[1e-7 1e4],'multiinicios',starts,'evaluaciones',evaluations, ...
    'flags',flags,'log_z_alpha',chosen);
    function [value,energy]=metrics(x)
        z0=exp(x(1)); alpha0=adm*exp(x(2));
        H=(D*D+z0*I)\D; Y=alpha0*H;
        u=(I+Z*(c.Yr+Y))\c.er; q=Y*u; is=(c.Yr+Y)*u;
        vc=u-D*q/alpha0;
        energy=dot(q,q)/(2*c.omega*alpha0)+alpha0*dot(vc,vc)/(2*c.omega*z0);
        value=dot(is,is)/baseline; evaluations=evaluations+1;
    end
    function value=objective(x), [value,~]=metrics(x); end
    function [ineq,eq]=constraint(x)
        [~,energy]=metrics(x); ineq=energy/cap-(1-1e-8); eq=[];
    end
    function consider(x,value)
        [~,energy]=metrics(x);
        if value<best && energy<=cap*(1+1e-9)
            best=value; chosen=x;
        end
    end
end
