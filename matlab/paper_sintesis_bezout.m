function r = paper_sintesis_bezout(c)
%PAPER_SINTESIS_BEZOUT Anulador temporal -> certificado -> Cauer por divisiones.
% Foster usa raices del denominador SOLO para extraer la realizacion alternativa.
N=c.N; A=-c.D*c.D; m=numel(c.n); E=zeros(N,m); x=c.u;
for k=1:m, E(:,k)=x; x=A*x; end
scale=vecnorm(E); p=[1;flipud((E./scale)\(-x)./scale')];
W=c.D*E; scaleW=vecnorm(W); h=flipud((W./scaleW)\c.qstar./scaleW');
p=p'; h=h'; den=conv([1 0],polyder(p))/m;
[~,num]=deconv(conv(den,h),p); num=trim(num);
H0=-bezout(den,p); B=bezout(den,num); L=chol(H0,'lower');
M=L\B/L'; M=(M+M')/2;
t=max(0,max(diag(M)+sum(abs(M),2)-abs(diag(M))));
% Interior estricto para no atribuir optimalidad de componentes a este certificado.
t=t+.01*max([1,t,norm(c.i)/norm(c.u)]);
num=[zeros(1,numel(p)-numel(num)),num]+t*p;
ascN=fliplr(num); ascD=fliplr(den);
yn=zeros(1,2*numel(ascN)); yd=zeros(1,2*numel(ascD)-1);
yn(2:2:end)=ascN.*(-1).^(0:numel(ascN)-1);
yd(1:2:end)=ascD.*(-1).^(0:numel(ascD)-1);
while yn(1)==0 && yd(1)==0, yn(1)=[]; yd(1)=[]; end
a=fliplr(yn); b=fliplr(yd); ladder=[];
for k=1:2*m+2
    a=trim(a); b=trim(b); assert(numel(a)==numel(b)+1,'Paper:CauerGrado','Grado incompatible.');
    value=a(1)/b(1); assert(value>0,'Paper:CauerSigno','Elemento no positivo.');
    ladder(end+1)=value/c.omega; %#ok<AGROW>
    rem=a-value*[b 0];
    if norm(rem,inf)<1e-9*max(norm(a,inf),1), break; end
    a=b; b=trim(rem);
end
assert(numel(ladder)<=2*m,'Paper:CauerOrden','Extraccion no terminada.');
mu=sort(real(roots(den))); alpha=-polyval(num,mu)./polyval(polyder(den),mu);
assert(all(alpha>-1e-8),'Paper:FosterSigno','Residuo negativo.');
branches=[0,t/c.omega,0];
for k=1:numel(mu)
    if alpha(k)<1e-12, continue; end
    if abs(mu(k))<1e-8, branches(end+1,:)=[1/(c.omega*alpha(k)),0,0]; %#ok<AGROW>
    else, branches(end+1,:)=[1/(c.omega*alpha(k)),alpha(k)/(c.omega*mu(k)),0]; end %#ok<AGROW>
end
r=struct('cauer',struct('ladder',ladder),'foster',struct('ramas',branches), ...
    'p',p,'h',h,'den',den,'num',num,'t',t, ...
    'certificado_max_eig',max(eig(B-t*H0)), ...
    'error_cierre',norm(E*flipud(p(2:end)')+x)/max(norm(x),1));
assert(r.certificado_max_eig<1e-7 && r.error_cierre<1e-8);
end
function z=trim(z)
while numel(z)>1 && abs(z(1))<1e-11*max(norm(z,inf),realmin), z(1)=[]; end
end
function B=bezout(a,b)
n=max(numel(a),numel(b))-1; a=fliplr(a); b=fliplr(b);
a(end+1:n+1)=0; b(end+1:n+1)=0; B=zeros(n);
for i=1:n
    for j=0:i-1
        value=a(i+1)*b(j+1)-b(i+1)*a(j+1);
        for k=0:i-j-1, B(i-k,j+k+1)=B(i-k,j+k+1)+value; end
    end
end
end
