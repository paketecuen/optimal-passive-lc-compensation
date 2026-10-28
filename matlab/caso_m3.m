function c = caso_m3(N)
%CASO_M3 Carga pasiva del certificado exacto; tiempo tau=2*pi*50*t.
arguments, N (1,1) double = 31, end
tau=(0:N-1)'*2*pi/N; D=derivador_periodico(N);
u=sqrt(2)*(4*cos(tau).^3+2*cos(tau).^2-2*cos(tau)-1);
u=u-mean(u); v=D*u;
aug=[D,ones(N,1);ones(1,N),0]; temp=aug\[u;0]; ju=temp(1:N);
G=.8; c0=11/20; gamma0=17/10; alpha=3/4; z=6;
i=G*u+c0*v+gamma0*ju+alpha*((D*D+z*eye(N))\v);
c=struct('u',u,'i',i,'D',D,'tau',tau,'omega',2*pi*50,'G',G, ...
    'C_normalizada',c0,'gamma_L',gamma0,'alpha_serie',alpha,'z_serie',z);
end
