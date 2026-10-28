function fit = ajustar_banco_medido(U,I,f,z,metodo)
%AJUSTAR_BANCO_MEDIDO Amplitudes positivas con una sintonia FISICA fijada.
% El filtro es temporal (I+lambda A^2)^-1; no usa FFT ni modos armonicos.
% No aplica una formula m3 exacta a medidas cuyo cierre es solo aproximado.
arguments
    U (:,:) double
    I (:,:) double
    f (1,1) double
    z (1,1) double
    metodo (1,1) string {mustBeMember(metodo,["un_periodo","promedio","regularizado"])}
end
N=size(U,1); P=size(U,2); D=derivador_periodico(N); A=-D*D;
if metodo=="un_periodo", u=U(:,1); i=I(:,1); else, u=mean(U,2); i=mean(I,2); end
u=u-mean(u); i=i-mean(i); raw=u; lambda=0; noise=sqrt(mean(var(U,0,2)))/sqrt(P);
if metodo=="regularizado" && noise>1e-10*sqrt(mean(u.^2))
    penalty=A*A; target=.9*noise; low=0; high=1e-12;
    while norm(u-(eye(N)+high*penalty)\u)/sqrt(N)<target && high<1
        high=high*10;
    end
    for k=1:35
        middle=(low+high)/2; candidate=(eye(N)+middle*penalty)\raw;
        if norm(raw-candidate)/sqrt(N)<target, low=middle; else, high=middle; end
    end
    lambda=high; u=(eye(N)+lambda*penalty)\raw;
end
ratio=f/50; aug=[D,ones(N,1);ones(1,N),0]; ju=aug\[u;0]; ju=ju(1:N);
w=(D*D+(z/ratio^2)*eye(N))\(D*u);
Phi=[ratio*(D*u),ju/ratio,w/ratio];
cone=cono_temporal(Phi,i);
fit=struct('coef',cone.alpha,'z',z,'metodo',metodo,'u_ajustada',u, ...
    'u_promedio',raw,'i_ajustada',i,'lambda',lambda,'sigma_promedio',noise, ...
    'correccion_u_rms',norm(u-raw)/sqrt(N),'I_predicha',sqrt(cone.I2), ...
    'kkt',cone.kkt_relativo,'q_predicha',cone.q);
end
