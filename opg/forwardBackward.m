function [Po] = forwardBackward(v,p,a,b)

% v - sekvenca opservacija
% p - pocetne verovatnoce
% a - verovatnoce prelaska iz stanja u stanje
% b - verovatnoce opservacije

% Po - verovatnoca

N = length(p);
T = length(v);

delta(1,:) = p.*b(:,v(1))';
psi(1,:) = zeros(1,N);
alpha = zeros(T,N);

%Inicijalizacija
alpha(1,:) = p.*b(:,v(1))';

%Indukcija
for t=2:T
    for j=1:N
        alpha(t,j) = sum(alpha(t-1,:).*a(:,j)');
        alpha(t,j) = alpha(t,j)*b(j,v(t));
    end
end

% Terminacija
Po = sum(alpha(T,:));
end