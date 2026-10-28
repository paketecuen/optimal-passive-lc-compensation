function f = familia_optima_m3(u,i,D,omega)
%FAMILIA_OPTIMA_M3 Bancos C || L || (LC serie) que realizan la proyeccion Q.
% Marco temporal H={v=Du,j=Ju,Av}; sus coordenadas son cocientes de
% trivectores. No se trata de una sola LC ni de una busqueda espectral.
arguments
    u (:,1) double
    i (:,1) double
    D (:,:) double
    omega (1,1) double {mustBePositive} = 2*pi*50
end
m=modelo_temporal_m3(u,D); q=-m.Q*(m.Q'*i/m.N);
H=[m.Du,m.Ju,m.A*m.Du]; M=m.Q'*H/m.N;
assert(rcond(M)>1e-10,'LCfamilia:Marco','Marco temporal degenerado.');
theta=M\(m.Q'*q/m.N); eta=theta(1); beta=theta(2); k=theta(3);
assert(abs(k)>1e-12*max(norm(theta),1),'LCfamilia:Paralelo', ...
    'La proyeccion ya pertenece al plano paralelo; esta carta supone k distinto de cero.');
% p(z)=z^3+p(3)z^2+p(2)z+p(1), en potencias ascendentes.
% W(z)=p(1)j+z Av+(z^2+p(3)z)v; w_z=W(z)/p(z).
% c=eta-k*(z+p(3)); gamma=beta-k*p(1)/z; alpha=k*p(z)/z.
pc=[eta-k*m.p(3),-k]; pg=[-k*m.p(1),beta]; pa=k*m.p;
cuts=[0;m.polos];
for poly={pc,pg,pa}
    r=roots(fliplr(poly{1})); r=real(r(abs(imag(r))<1e-9 & real(r)>0));
    cuts=[cuts;r]; %#ok<AGROW>
end
cuts=sort(cuts); cuts=cuts([true;diff(cuts)>1e-8*max(1,max(cuts))]);
cuts=[cuts;Inf]; rows=zeros(0,4);
for n=1:numel(cuts)-1
    lo=cuts(n); hi=cuts(n+1); probe=(lo+hi)/2;
    if isinf(hi), probe=lo+max(lo,1); end
    coeff=[polyval(fliplr(pc),probe);polyval(fliplr(pg),probe)/probe; ...
        polyval(fliplr(pa),probe)/probe];
    if all(coeff>0)
        rows(end+1,:)=[lo,hi,endpoint(lo),endpoint(hi)]; %#ok<AGROW>
    end
end
f=struct('modelo',m,'q_ideal',q,'marco',H,'coordenadas',theta, ...
    'omega',omega,'c_numerador',pc,'gamma_numerador',pg, ...
    'alpha_numerador',pa,'error_marco_A',norm(H*theta-q)/sqrt(m.N), ...
    'I_ideal_A',norm(i+q)/sqrt(m.N),'I_original_A',norm(i)/sqrt(m.N), ...
    'intervalos',array2table(rows,'VariableNames', ...
    {'z_min','z_max','Incluye_min','Incluye_max'}));

    function yes=endpoint(z)
        yes=false;
        if ~isfinite(z) || z<=0 || min(abs(z-m.polos))<1e-8, return; end
        a=[polyval(fliplr(pc),z);polyval(fliplr(pg),z)/z;polyval(fliplr(pa),z)/z];
        yes=all(a(1:2)>=-1e-9*max(norm(a),1)) && a(3)>0;
    end
end
