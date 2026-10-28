function D = derivador_periodico(N)
%DERIVADOR_PERIODICO Derivada en tau del interpolante periodico de N nodos.
% N impar, tau_j=2*pi*j/N. El modelo finito se declara, no se infiere de
% cualquier senal analogica que pase por los nodos. No se ejecuta FFT.
arguments
    N (1,1) double {mustBeInteger,mustBePositive}
end
assert(mod(N,2)==1 && N>=3,'LCm2:Nodos','Se requiere N impar, N>=3.');
d = (0:N-1)'-(0:N-1);
D = zeros(N);
mask = d~=0;
D(mask) = 0.5*(-1).^d(mask)./sin(pi*d(mask)/N);
end
