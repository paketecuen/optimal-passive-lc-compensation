function out = interpolar_periodico(values,tau)
%INTERPOLAR_PERIODICO Formula cardinal temporal del interpolante impar.
N=size(values,1); theta=tau(:)-(0:N-1)*2*pi/N;
den=sin(theta/2); K=sin(N*theta/2)./(N*den);
K(abs(den)<1e-13)=1; out=K*values;
end
