%% Doma?i zadatak iz prepoznavanja oblika - zadatak 3

clear;
close all;
clc;

%% Generisanje dve linearno separabilne klase dvodimenzinalnih oblika

N=500;

%Klasa 1

M1 = [0; 4]; 
S1 = [4 1.1;1.1 2];

[F1,L1]=eig(S1);
T1=F1*L1^(1/2);

pomv = rand(1,N);
X1 = (T1*randn(2,N)+M1*ones(1,N));

figure(8);
plot(X1(1,:),X1(2,:),'ro');
hold all

%Klasa 2

M2 = [6; 0];
S2 = [3 0.8;0.8 0.5];

[F2,L2]=eig(S2);
T2=F2*L2^(1/2);

pomv = rand(1,N);
X2 = (T2*randn(2,N)+M2*ones(1,N));

figure(8);
plot(X2(1,:),X2(2,:),'bo');
xlabel('X1');
ylabel('X2');
legend('K1','K2');
title('Odbirci klasa');

%% 

%podela skupova na trenirajuci i testirajuci skup

X1_tr=X1(:,1:round(0.7*N));
X1_test = X1(:,round(0.7*N)+1:N);

X2_tr=X2(:,1:round(0.7*N));
X2_test = X2(:,round(0.7*N)+1:N);

%procena srednje vrednosti i kovar. matrica na osnovu obucavajuceg skupa

M1_est=mean(X1_tr');
M2_est=mean(X2_tr');

S1_est=cov(X1_tr');
S2_est=cov(X2_tr');
%% Linearni klasifikator - hold-out metoda

errornum_s_min=0;
s_opt=0;
v0_s_opt=0;
for s=0:0.01:1
    V=(s*S1_est+(1-s)*S2_est)^(-1)*(M2_est'-M1_est');
    y1=zeros(1,length(X1_test));
    y2=zeros(1,length(X2_test));
    for i=1:length(X1_test)
        y1(1,i) = V'*X1_test(:,i);
        y2(1,i) = V'*X2_test(:,i);
    end
    v0min=-max(max(y1),max(y2));
    v0max=-min(min(y1),min(y2));
    delta_v0=0.01;
    errornum_min=0;
    v0_opt=0;
    for v0=v0min:delta_v0:v0max
        errornum=0;
        for i=1:length(X1_test)
            if y1(1,i)>-v0
                errornum=errornum+1;
            end
            if y2(1,i)<-v0
                errornum=errornum+1;
            end
        end  
        if v0==v0min
            errornum_min=errornum;
            v0_opt=v0;
        else
            if errornum<errornum_min
                errornum_min=errornum;
                v0_opt=v0;
            end
        end
    end
    if s==0
        errornum_s_min=errornum_min;
        s_opt=s;
        v0_s_opt=v0_opt;
    else
        if errornum_min<errornum_s_min
            errornum_s_min=errornum_min;
            s_opt=s;
            v0_s_opt=v0_opt;   
        end
    end 
end

%ponovo racunanje parametara sa optimalnim s koje smo dobili iznad

V=(s_opt*S1_est+(1-s_opt)*S2_est)^(-1)*(M2_est'-M1_est');
y1=zeros(1,length(X1_test));
y2=zeros(1,length(X2_test));

for i=1:length(X2_test)
    y1(1,i) = V'*X1_test(:,i);
    y2(1,i) = V'*X2_test(:,i);
end

v0min=-max(max(y1),max(y2));
v0max=-min(min(y1),min(y2));
delta_v0=0.01;
errornum_min=0;
v0_opt=0;
for v0=v0min:delta_v0:v0max
    errornum=0;
    for i=1:length(X1_test)
        if y1(1,i)>-v0
            errornum=errornum+1;
        end
        if y2(1,i)<-v0
            errornum=errornum+1;
        end
    end  
    if v0==v0min
        errornum_min=errornum;
        v0_opt=v0;
    else
        if errornum<errornum_min
            errornum_min=errornum;
            v0_opt=v0;
        end
    end
end


figure(9);
plot(X1(1,:),X1(2,:),'ro')
hold on
plot(X2(1,:),X2(2,:),'bo')

x = -10:0.1:12;
y = -4:0.1:8;

for i = 1: length(x)
    for j = 1: length(y)
        X = [x(i);  y(j)];
        h(i,j)= V'*X+v0_opt;   
    end
end

contour(x,y,h',[0 0],'r--','linewidth',1.5); 
legend('K1','K2','klasifikaciona linija');
title('Linearni klasifikator');
xlim([-10 12]);
ylim([-4 8]);

%% Linearni klasifikator - metoda metoda 

Gama=ones(2*N,1);
U=[-1*ones(1,N) 1*ones(1,N);-X1 X2];
W=(U*U')^(-1)*U*Gama;
v0 = W(1);
v = W(2:3);

xpReal = [-10 12];
figure(10);
hold on
plot(X1(1,:),X1(2,:),'ro');
plot(X2(1,:),X2(2,:),'bo');
plot(xpReal,(-v0-v(1)*xpReal)/v(2),'y--','linewidth',1.5);

Gama = [1.5*ones(N,1); ones(N,1)];
U=[-1*ones(1,N) 1*ones(1,N);-X1 X2];
W=inv((U*U'))*U*Gama;
v0 = W(1);
v = W(2:3);

xpReal = [-10 12];
plot(xpReal,(-v0-v(1)*xpReal)/v(2),'m--','linewidth',1.5);

Gama = [ones(N,1); 1.5*ones(N,1)];
U=[-1*ones(1,N) 1*ones(1,N);-X1 X2];
W=inv((U*U'))*U*Gama;
v0 = W(1);
v = W(2:3);

xpReal = [-10 12];
plot(xpReal,(-v0-v(1)*xpReal)/v(2),'c--','linewidth',1.5);

hold off
title('Linearni klasifikator metodom zeljenog izlaza');
legend('K1','K2');
xlim([-10 12]);
ylim([-4 8]);

%% Kvadratni klasifikator metodom zeljenog izlaza za dve linearno neseparabilne klase

clear ;
close all;
clc;


%% Generisanje odbiraka 

N=500;

%Klasa 1

M1=[0;0];
R1=rand(1,N);
alfa1=2*pi*rand(1,N);
X1=[R1.*cos(alfa1);R1.*sin(alfa1)]+M1*ones(1,N);

figure(1);
plot(X1(1,:),X1(2,:),'ro');

%Klasa 2

M2=[0;0];
R_u=3;
R_d=2;
R2=R_d*rand(1,N)+R_u;
alfa2=2*pi*rand(1,N);
X2=[R2.*cos(alfa2);R2.*sin(alfa2)]+M2*ones(1,N);
hold on
plot(X2(1,:),X2(2,:),'bo');
legend('K1','K2')
hold off

%% Projektovanje klasifikatora

Gama=ones(2*N,1);
U=[-1*ones(1,N) 1*ones(1,N)  ;  -X1 X2  ;  -X1(1,:).^2 X2(1,:).^2  ;  -2*X1(1,:).*X1(2,:) 2*X2(1,:).*X2(2,:)  ;  -X1(2,:).^2 X2(2,:).^2];
W=(U*U')^(-1)*U*Gama;
v0=W(1);
v=W(2:3);
Q=[W(4) W(5);W(5) W(6)];

syms xp yp
[xp,yp,~,~]=solve(v0+xp*v(1)+yp*v(2)+xp^2*Q(1)+2*xp*yp*Q(2)+yp^2*Q(4),xp,yp,'returnConditions',true);
z=-4:0.001:4;
xp=eval(xp);
xp=[xp(1,:) fliplr(xp(2,:))];
xpReal=xp(imag(xp)==0);
yp=eval(yp);
yp=[yp(1,:) fliplr(yp(2,:))];
ypReal=yp(imag(xp)==0);

figure(11);
plot(X1(1,:),X1(2,:),'ro');
hold on
plot(X2(1,:),X2(2,:),'bo');
plot(xpReal,ypReal,'g--','LineWidth',1.5);
legend('K1','K2','Diskriminaciona linija')
hold off

