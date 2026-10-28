function r = paper_evaluar_red(c,b,opciones)
%PAPER_EVALUAR_RED Circuito periodico por operadores temporales en su cierre.
% Banco Foster: filas [L_H,C_F,R_ohm], cero=elemento ausente.
% O Cauer-I: b.ladder, valores SI, alterna C shunt/L serie.
arguments
    c (1,1) struct
    b (1,1) struct
    opciones.Trazas (1,1) logical = true
    opciones.Puntos (1,1) double = 4096
end
D=c.omega*c.Dr; I=eye(size(D)); Z=c.Rg*I+c.Lg*D;
if isfield(b,'ladder')
    tail=cell(numel(b.ladder),1);
    for k=numel(tail):-1:1
        tail{k}=b.ladder(k)*D;
        if k<numel(tail), tail{k}=tail{k}+tail{k+1}\I; end
    end
    startsC=~isfield(b,'ladderStart') || b.ladderStart=="C";
    Y=tail{1}; if ~startsC, Y=Y\I; end
    elements=numel(tail);
else
    Y=zeros(size(D)); branches=cell(size(b.ramas,1),1); elements=0;
    for k=1:size(b.ramas,1)
        L=b.ramas(k,1); C=b.ramas(k,2); R=b.ramas(k,3);
        if L==0 && C==0, branches{k}=zeros(size(D)); continue; end
        if C==0
            branches{k}=(L*D+R*I)\I;
        elseif L==0
            branches{k}=(I+R*C*D)\(C*D);
        else
            branches{k}=(I+R*C*D+L*C*D*D)\(C*D);
        end
        Y=Y+branches{k}; elements=elements+(L>0)+(C>0);
    end
end
u=(I+Z*(c.Yr+Y))\c.er; il=c.Yr*u; q=Y*u; is=il+q;
r=struct('Is_RMS_A',norm(is),'U_RMS_V',norm(u), ...
    'PF_fuente',dot(c.er,is)/(norm(c.er)*norm(is)), ...
    'PF_bornes',dot(u,is)/(norm(u)*norm(is)), ...
    'P_carga_W',dot(u,il),'P_banco_W',dot(u,q), ...
    'P_red_W',c.Rg*dot(is,is),'Elementos',elements);
if ~opciones.Trazas, return; end
tau=(0:opciones.Puntos-1)'*2*pi/opciones.Puntos;
K=interpolar_periodico(c.V,tau); vv=[]; ii=[]; energy=zeros(size(tau));
if isfield(b,'ladder')
    voltage=u; current=q;
    for k=1:numel(b.ladder)
        if (mod(k,2)==1)==startsC
            vi=voltage; ji=b.ladder(k)*D*voltage;
            energy=energy+.5*b.ladder(k)*(K*vi).^2; current=current-ji;
        else
            ji=current; vi=b.ladder(k)*D*current;
            energy=energy+.5*b.ladder(k)*(K*ji).^2; voltage=voltage-vi;
        end
        vv=[vv,K*vi]; ii=[ii,K*ji]; %#ok<AGROW>
    end
else
    for k=1:size(b.ramas,1)
        L=b.ramas(k,1); C=b.ramas(k,2); R=b.ramas(k,3);
        ji=branches{k}*u; vl=L*D*ji; vc=zeros(size(u));
        if C>0, vc=u-vl-R*ji; end
        energy=energy+.5*L*(K*ji).^2+.5*C*(K*vc).^2;
        vv=[vv,K*vl,K*vc]; ii=[ii,K*ji]; %#ok<AGROW>
    end
end
r.E_media_J=mean(energy); r.E_max_J=max(energy);
r.V_elemento_pico_V=max([0;abs(vv(:))]); r.I_rama_pico_A=max([0;abs(ii(:))]);
r.traza=struct('tau',tau,'u',K*u,'is',K*is,'il',K*il,'q',K*q, ...
    'v_elementos',vv,'i_ramas',ii,'energia',energy);
r.q_nodos=c.V*q; r.u_nodos=c.V*u; r.is_nodos=c.V*is;
end
