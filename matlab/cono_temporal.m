function r = cono_temporal(Phi,i)
%CONO_TEMPORAL Proyeccion al cono positivo por todas sus caras independientes.
% Las coordenadas se obtienen del marco reciproco (QR escalada), no por
% recorte de la proyeccion libre. No requiere Optimization Toolbox.
N=numel(i); k=size(Phi,2); assert(k<=10,'LC:Cono','Demasiadas ramas para enumerar caras.');
best=mean(i.^2); alpha=zeros(k,1); active=[]; scale=sqrt(mean(Phi.^2,1));
for mask=1:2^k-1
    S=find(bitget(mask,1:k));
    if any(scale(S)<realmin), continue; end
    V=Phi(:,S)./scale(S); [U,R]=qr(V,0);
    if size(R,1)<numel(S) || rcond(R)<1e-10, continue; end
    c=-(R\(U'*i)); tol=1e-10*max(norm(i),realmin);
    if any(c < -tol), continue; end
    c=max(c,0); cand=i+V*c; value=mean(cand.^2);
    if value < best-32*eps*max(mean(i.^2),realmin)
        best=value; alpha=zeros(k,1); alpha(S)=c./scale(S)'; active=S;
    end
end
q=Phi*alpha; residual=i+q; grad=Phi'*residual/N;
den=max(scale',realmin)*max(sqrt(mean(i.^2)),realmin);
kkt=max([0;-grad./den;abs(grad(alpha>0))./den(alpha>0)]);
r=struct('alpha',alpha,'q',q,'residuo',residual,'I2',mean(residual.^2), ...
    'activos',active,'kkt_relativo',kkt);
end
