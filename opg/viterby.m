function [x,Po] = viterby(v,p,a,b)

% v - sekvenca opservacija
% p - pocetne verovatnoce
% a - verovatnoce prelaska iz stanja u stanje
% b - verovatnoce opservacije

% Po - verovatnoca
% x - sekvenca stanja

N = length(p);
T = length(v);

% Inicijalizacija
delta(1,:) = p.*b(:,v(1))';
psi(1,:) = zeros(1,N);

% Rekurzija
for t=2:T
    for j=1:N
        [delta(t,j), psi(t,j)] = max(delta(t-1,:).*a(:,j)');
        delta(t,j)=delta(t,j)*b(j,v(t));
    end
end

% Terminacija
[Po(t),x(T)] = max(delta(T,:));

% Razmotavanje
for t=T-1:-1:1
    x(t)=psi(t+1, x(t+1));
    Po(t) = delta(t,x(t));
end

end