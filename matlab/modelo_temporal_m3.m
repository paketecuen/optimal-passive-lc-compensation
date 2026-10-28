function m = modelo_temporal_m3(u,D,opciones)
%MODELO_TEMPORAL_M3 Marco RMS de Q=span(Du,A Du,A^2 Du), A=-D^2.
% Se genera el mismo espacio con A^-1 para evitar derivadas altas.
% Ortonormalizacion temporal, sin FFT ni diagonalizacion modal.
arguments
    u (:,1) double {mustBeReal,mustBeFinite}
    D (:,:) double {mustBeReal,mustBeFinite}
    opciones.ToleranciaCierre (1,1) double {mustBePositive} = 1e-9
end
N=numel(u); rmsv=@(x) norm(x)/sqrt(N);
assert(isequal(size(D),[N,N]),'LCm3:Dimensiones','D incompatible.');
assert(rmsv(u)>0 && abs(mean(u))<1e-11*rmsv(u),'LCm3:DC','Se requiere tension AC no nula.');
assert(norm(D+D','fro')<1e-11*max(norm(D,'fro'),1), ...
    'LCm3:Derivador','D debe ser antisimetrico.');
assert(norm(D*ones(N,1))<1e-11*max(norm(D,'fro'),1), ...
    'LCm3:Derivador','D debe anular constantes.');
A=-D*D; v=D*u; Q=zeros(N,3); x=v;
augA=[A,ones(N,1);ones(1,N),0];
for k=1:3
    for twice=1:2, x=x-Q(:,1:k-1)*(Q(:,1:k-1)'*x/N); end
    assert(rmsv(x)>1e-11*max(rmsv(v),1),'LCm3:Rango','Dimension menor que tres o mal resuelta.');
    Q(:,k)=x/rmsv(x);
    temp=augA\[Q(:,k);0]; x=temp(1:N);
end
T=Q'*A*Q/N; T=(T+T')/2;
closure=norm(A*Q-Q*T,'fro')/max(norm(A*Q,'fro'),realmin);
assert(closure<=opciones.ToleranciaCierre,'LCm3:Cierre', ...
    'No hay cierre m3: no sustituirlo por un ajuste sin declarar.');
p=[-det(T),(trace(T)^2-trace(T*T))/2,-trace(T),1];
vc=Q'*v/N;
W=Q*[(T*T+p(3)*T+p(2)*eye(3))*vc,(T+p(3)*eye(3))*vc,vc];
poles=sort(roots(fliplr(p)));
assert(all(abs(imag(poles))<1e-8) && all(real(poles)>0), ...
    'LCm3:Polos','El cierre debe tener tres raices positivas.');
poles=real(poles);
assert(min(diff(poles))>1e-8*max(poles),'LCm3:Polos','Raices no resueltas.');
aug=[D,ones(N,1);ones(1,N),0]; sol=aug\[u;0]; ju=sol(1:N);
assert(rmsv(D*ju-u)<1e-9*rmsv(u),'LCm3:Primitiva','Primitiva incompatible.');
m=struct('u',u,'D',D,'A',A,'Du',v,'Ju',ju,'Q',Q,'T',T, ...
    'p',p,'W',W,'polos',poles,'cierre_relativo',closure,'N',N);
end
