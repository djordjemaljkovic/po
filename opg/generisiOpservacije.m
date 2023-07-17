function [s, v] = generisiOpservacije(N,p,a,b)

% N = duzina sekvence
% p - pocetne verovatnoce
% a - verovatnoca prelaska iz stanja u stanje
% b - verovatnoce opservacija

% s - sekvenca stanja
% v - sekvenca opservacija
% 1 = bela kuglica
% 2 = crna kuglica
% 3 = zelena kuglica

s = zeros(1,N+1);       %vektor stanja
v = zeros(1,N);         %vektor opservacija

P = cumsum(p);
A = cumsum(a,2);        %matrica tranzicija
B = cumsum(b,2);        %matrica opservacija


temp = rand();
s(1) = find(temp<P,1);

for i=1:N
    % Generisanje opservacije
    temp = rand();
    v(i) = find(temp<B(s(i),:),1);
    
    % Generisanje sledeceg stanja
    temp = rand();
    s(i+1) = find(temp<A(s(i),:),1);
end
s = s(1:end-1);
end