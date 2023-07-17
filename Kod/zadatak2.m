%% Doma?i zadatak iz prepoznavanja oblika - zadatak 2

clear; 
close all;
clc;

%% Generisanje odbiraka i prikaz dijagrama

% Klasa 1 : 0.6*N11 + 0.4*N12
% Klasa 2 : 0.45*N21 + 0.55*N22

N=500;

%Klasa 1
P11 = 0.6;
P12 = 0.4;
M11 = [1; 1];
S11 = [4 1.1;1.1 2];
M12 = [6; 4]; 
S12 = [3 -0.8;-0.8 1.5];

[F11,L11]=eig(S11);
T11=F11*L11^(1/2);
[F12,L12]=eig(S12);
T12=F12*L12^(1/2);

pomv = rand(1,N);
X1 = (T11*randn(2,N) + M11*ones(1,N)).*[pomv<=P11;pomv<=P11]+(T12*randn(2,N)+M12*ones(1,N)).*[pomv>P11;pomv>P11];



%Klasa 2
P21 = 0.55;
P22 = 0.45;
M21 = [7; -4]; S21 = [2 1.1; 1.1 4];
M22 = [6; 0]; S22 = [3 0.8;0.8 0.5];

[F21,L21]=eig(S21);
T21=F21*L21^(1/2);
[F22,L22]=eig(S22);
T22=F22*L22^(1/2);

pomv = rand(1,N);
X2 = (T21*randn(2,N)+M21*ones(1,N)).*[pomv<=P21;pomv<=P21]+(T22*randn(2,N)+M22*ones(1,N)).*[pomv>P21;pomv>P21];




%% Funkcija gustine verovatnoce i diskriminaciona funkcija

x = -6:0.1:12;
y = -10:0.1:8;
for i = 1: length(x)
    for j = 1: length(y)
        X = [x(i);  y(j)];
        f11 = 1/(2*pi*det(S11)^0.5)*exp(-0.5*(X-M11)'*inv(S11)*(X-M11));
        f12 = 1/(2*pi*det(S12)^0.5)*exp(-0.5*(X-M12)'*inv(S12)*(X-M12));
        f21 = 1/(2*pi*det(S21)^0.5)*exp(-0.5*(X-M21)'*inv(S21)*(X-M21));
        f22 = 1/(2*pi*det(S22)^0.5)*exp(-0.5*(X-M22)'*inv(S22)*(X-M22));
        f1(i,j) = P11*f11 + P12*f12;
        f2(i,j) = P21*f21 + P22*f22;
        h(i,j) = -log(f1(i,j)) + log(f2(i,j));
         
    end
end

figure(2);
mesh(x,y,f1');
hold on;
mesh(x,y,f2'); 
hold off;
title('Funkcije gustine verovatnoce');
xlim([-6 12]);
ylim([-10 8]);

%% d2 krive

% Klasa 1
d2 = [1 4 9];
prag = max(max(f1))*exp(-0.5*d2);
figure(1)
plot(X1(1,:),X1(2,:),'rh');
hold all
contour(x,y,f1',prag,'m');
hold all

% Klasa 2
prag = max(max(f2))*exp(-0.5*d2);
figure(1)
plot(X2(1,:),X2(2,:),'bh');
xlabel('X1');
ylabel('X2');
contour(x,y,f2',prag,'b');
title('d^2 krive u prostoru oblika')
legend('K1','d2=1, 4, 9','K2', 'd2=1, 4, 9');
xlim([-6 12]);
ylim([-10 8]);

%% Bajesov klasifikator minimalne greske

%uzimamo vrednost za h iz prethodnog dela zadatka 2

%P1=P2=0.5, pretpostavka

figure(1);
plot(X1(1,:),X1(2,:),'rh');
hold all;
contour(x,y,f1',prag,'m');
plot(X2(1,:),X2(2,:),'bh');
xlabel('X1');
ylabel('X2');
contour(x,y,f2',prag,'b');
contour(x,y,h',[0 0],'g--','linewidth',1.5); 
legend('K1','K2','','','klasifikaciona linija');
title('Bajesov klasifikator minimalne greske')
xlim([-6 12]);
ylim([-10 8]);

% verovatnoca greske Bajesovog klasifikatora

Eps1 = 0;                       %verovatnoca greske prvog tipa
Eps2 = 0;                       %verovatnoca greske drugog tipa
for i = 1:length(x)
    for j = 1:length(y)
        if h(i,j) < 0
            Eps2 = Eps2 + f2(i,j)*0.1*0.1; 
        else 
            Eps1 = Eps1 + f1(i,j)*0.1*0.1;
        end
    end
end

Eps1  

Eps2

Eps=0.5*(Eps1+Eps2)



%% Bajesov test minimalne cene


c11 = 0;            % unapred usvojene cene
c22 = 0;
c12 = 0.9;
c21 = 0.1;
pragnovi = -log((c12-c22)/(c21-c11));
figure(2);
hold all;
plot(X1(1,:),X1(2,:),'rh');
plot(X2(1,:),X2(2,:),'bh');
contour(x,y,f1',prag,'m');
contour(x,y,f2',prag,'b');
contour(x,y,h',[0 0],'g--','linewidth',1.5); 
contour(x,y,h',[pragnovi pragnovi],'k--','linewidth',1.5);
xlabel('X1');
ylabel('X2');
legend('K1','K2','','','minimalna greska','minimalna cena');
title('Bajesov klasifikator minimalne cene')
xlim([-6 12]);
ylim([-10 8]);

% verovatnoca greske Bajesovog klasifikatora

Eps1 = 0;                       %verovatnoca greske prvog tipa
Eps2 = 0;                       %verovatnoca greske drugog tipa
for i = 1:length(x)
    for j = 1:length(y)
        if h(i,j) < pragnovi
            Eps2 = Eps2 + f2(i,j)*0.1*0.1; 
        else 
            Eps1 = Eps1 + f1(i,j)*0.1*0.1;
        end
    end
end

Eps1

Eps2

%% Wald-ov sekvencijalni test

eps1=1e-5;          %unapred usvojene greske
eps2=1e-5;
a=-log((1-eps1)/eps2);
b = -log(eps1/(1-eps2));

figure(3);
hold all

%Prva klasa
for i = 1:N
    kraj=0;
    Sm = 0;
    SM = [0];
    ind=randperm(N);
    K1(1,:)=X1(1,ind);
    K1(2,:)=X1(2,ind);
    j=0;
    while kraj==0
        j=j+1;
        odb=[K1(1,j);K1(2,j)];
        f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
        f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
        f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
        f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
        f1= P11*f11 + P12*f12;
        f2= P21*f21 + P22*f22;
        hx=-log(f1)+ log(f2);
        Sm = Sm+hx;
        SM=[SM Sm];
        if Sm<a || Sm>b
            kraj=1;
            t=0:1:(length(SM)-1);
            plot(t,SM,'r');
        end 
    end
end

%Druga klasa
for i = 1:N
    kraj=0;
    Sm = 0;
    SM = [0];
    ind=randperm(N);
    K2(1,:)=X2(1,ind);
    K2(2,:)=X2(2,ind);
    j=0;    
    while kraj==0  
        j=j+1;
        odb=[K2(1,j);K2(2,j)];
        f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
        f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
        f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
        f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
        f1= P11*f11 + P12*f12;
        f2= P21*f21 + P22*f22;
        hx=-log(f1)+ log(f2);
        Sm = Sm+hx;
        SM=[SM Sm];  
        if Sm<a || Sm>b
            kraj=1;
            t=0:1:(length(SM)-1);
            plot(t,SM,'b');
        end    
    end
end

plot([0,6],[a,a],'g--','LineWidth',1.5)
plot([0,6],[b,b],'g--','LineWidth',1.5)
title('Wald-ov sekvencijalni test');
grid on
hold off

%% Prva klasa - zavisnost od verovatnoce greske prvog tipa

e1=logspace(-10,0,1000);
m1=[];

for k=1:length(e1)-1
    a=-log((1-e1(k))/eps2);
    b = -log(e1(k)/(1-eps2));
    brOdb=0;
    
    %Prva klasa
    for i = 1:N
        kraj=0;
        Sm = 0;
        SM = [0];
        ind=randperm(N);
        K1(1,:)=X1(1,ind);
        K1(2,:)=X1(2,ind);
        j=0;
        while kraj==0
            j=j+1;
            odb=[K1(1,j);K1(2,j)];
            f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
            f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
            f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
            f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
            f1= P11*f11 + P12*f12;
            f2= P21*f21 + P22*f22;
            hx=-log(f1)+ log(f2);
            Sm = Sm+hx;
            SM=[SM Sm];
            if Sm<a || Sm>b
                kraj=1;
                brOdb=brOdb+j;
            end
        end
    end
    prosek=brOdb/N;
    m1=[m1 prosek]; 
end


figure(4);
semilogx(e1(1:length(e1)-1),m1);
xlabel('Greska prvog tipa'); 
title('Zavisnost broja potrebnih odbiraka prve klase od usvojene verovatnoce gresaka prvog tipa');


%% Prva klasa - zavisnost od verovatnoce greske drugog tipa

e2=logspace(-10,0,1000);
m2=[];

for k=1:length(e2)-1
    a=-log((1-eps1)/e2(k));
    b = -log(eps1/(1-e2(k)));
    brOdb=0;

    %Prva klasa
    for i = 1:N
        kraj=0;
        Sm = 0;
        SM = [0];
        ind=randperm(N);
        K1(1,:)=X1(1,ind);
        K1(2,:)=X1(2,ind);
        j=0;
        while kraj==0
            j=j+1;
            odb=[K1(1,j);K1(2,j)];
            f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
            f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
            f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
            f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
            f1= P11*f11 + P12*f12;
            f2= P21*f21 + P22*f22;
            hx=-log(f1)+ log(f2);
            Sm = Sm+hx;
            SM=[SM Sm];
            if Sm<a || Sm>b
                kraj=1;
                t=0:1:(length(SM)-1);
                brOdb=brOdb+j;
            end
        end
    end
    prosek=brOdb/N;
    m2=[m2 prosek]; 
end

figure(5);
semilogx(e2(1:length(e2)-1),m2);
xlabel('Greska drugog tipa');
title('Zavisnost broja potrebnih odbiraka prve klase od usvojene verovatnoce gresaka drugog tipa');

%% Druga klasa - zavisnost od verovatnoce greske prvog tipa

e1=logspace(-10,0,1000);
m1=[];

for k=1:length(e1)-1  
    a=-log((1-e1(k))/eps2);
    b = -log(e1(k)/(1-eps2));
    brOdb=0;

    %Druga klasa
    for i = 1:N
        kraj=0;
        Sm = 0;
        SM = [0];
        ind=randperm(N);
        K2(1,:)=X2(1,ind);
        K2(2,:)=X2(2,ind);
        j=0;
        while kraj==0
            j=j+1;
            odb=[K2(1,j);K2(2,j)];
            f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
            f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
            f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
            f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
            f1= P11*f11 + P12*f12;
            f2= P21*f21 + P22*f22;
            hx=-log(f1)+ log(f2);
            Sm = Sm+hx;
            SM=[SM Sm];
            if Sm<a || Sm>b
                kraj=1;
                brOdb=brOdb+j;
            end
        end
    end
    prosek=brOdb/N;
    m1=[m1 prosek];
end

figure(6);
semilogx(e1(1:length(e1)-1),m1);
xlabel('Greska prvog tipa');
title('Zavisnost broja potrebnih odbiraka druge klase od usvojene verovatnoce gresaka prvog tipa');

%% Druga klasa - zavisnost od verovatnoce greske drugog tipa

e2=logspace(-10,0,1000);
m2=[];

for k=1:length(e2)-1    
    a=-log((1-eps1)/e2(k));
    b = -log(eps1/(1-e2(k)));
    brOdb=0;

    %Prva klasa
    for i = 1:N
        kraj=0;
        Sm = 0;
        SM = [0];
        ind=randperm(N);
        K2(1,:)=X2(1,ind);
        K2(2,:)=X2(2,ind);
        j=0;
        while kraj==0
            j=j+1;
            odb=[K2(1,j);K2(2,j)];
            f11=P11/(2*pi*det(S11)^0.5)*exp(-0.5*(odb-M11)'*inv(S11)*(odb-M11));
            f12=P12/(2*pi*det(S12)^0.5)*exp(-0.5*(odb-M12)'*inv(S12)*(odb-M12));
            f21=P21/(2*pi*det(S21)^0.5)*exp(-0.5*(odb-M21)'*inv(S21)*(odb-M21));
            f22=P22/(2*pi*det(S22)^0.5)*exp(-0.5*(odb-M22)'*inv(S22)*(odb-M22));
            f1= P11*f11 + P12*f12;
            f2= P21*f21 + P22*f22;
            hx=-log(f1)+ log(f2);
            Sm = Sm+hx;
            SM=[SM Sm];
            if Sm<a || Sm>b
                kraj=1;
                t=0:1:(length(SM)-1);
                brOdb=brOdb+j;
            end
        end
    end
    prosek=brOdb/N;
    m2=[m2 prosek];
end

figure(7);
semilogx(e2(1:length(e2)-1),m2);
xlabel('Greska drugog tipa');
title('Zavisnost broja potrebnih odbiraka druge klase od usvojene verovatnoce gresaka drugog tipa');


















