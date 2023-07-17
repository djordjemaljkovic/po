function e=estimate(yy,lambda,tau,win,fs)

% yy - sekvenca
% lambda - parametar eksponencijalnog opadanja
% tau - blanking period
% win - sirina prozora
% fs - frekvencija odabiranja
% e - pitch frekvencija

j=1;
while yy(j)==0 && j<win
    j=j+1;
end

mt=yy(j); 
t=j; 
j=min(win,j+tau-1); 
while yy(j)<mt*exp(-lambda*(j-t)) && j<win
    j=j+1;
end

e=j-t;  
e = e/fs;
end