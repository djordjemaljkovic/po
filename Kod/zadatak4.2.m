%% Domaci zadatak iz prepoznavanja oblika - zadatak 4, zahtev 2

clear;
close all;
clc;

%% Generisanje 2 nelinearno separabilne klase

N=500;

%Klasa 1 

M1=[0.5;1.5];
R1=rand(1,N);
alfa1=pi*rand(1,N);
X1=[R1.*cos(alfa1);R1.*sin(alfa1)]+M1*ones(1,N);
figure(20);
plot(X1(1,:),X1(2,:),'o','color',[0.8500 0.3250 0.0980])

%Klasa 2 

M2=[0;0];
R_u=3;
R_d=2;
R2=R_d*rand(1,N)+R_u;
alfa2=pi*rand(1,N);
X2=[R2.*cos(alfa2);R2.*sin(alfa2)]+M2*ones(1,N);
hold on
plot(X2(1,:),X2(2,:),'o','color',[0.4660 0.6740 0.1880])
title('Odbirci klasa');
legend('K1','K2')
xlabel('X1');
ylabel('X2');
hold off

%% Metoda kvadratne dekompozicije


%% Nasumicna pocetna klasterizacija sa 2 klase 


pom=rand(2,N);

K1=[];
K2=[];

for i=1:N                   % skup X1
    if pom(1,i)<0.5
        K1=[K1 X1(:,i)];
    else
        K2=[K2 X1(:,i)];
    end
end

for i=1:N                   % skup X2
    if pom(2,i)<0.5
        K1=[K1 X2(:,i)];
    else
        K2=[K2 X2(:,i)];

    end
end

figure(21); 
plot(K1(1,:),K1(2,:),'o','color',[0.8500 0.3250 0.0980]);
hold on; 
plot(K2(1,:),K2(2,:),'o','color',[0.4660 0.6740 0.1880]);
title('Pocetna klasterizacija (nasumicno)');
xlabel('X1');
ylabel('X2');
legend('K1','K2');
hold off

%% Pocetna klasterizacija sa apriornim predznanjem sa 2 klase


pom=rand(2,N*3/10);

K1=X1(:,N*3/10+1:end);
K2=X2(:,N*3/10+1:end);

for i=1:N*3/10              % skup X1
    if pom(1,i)<0.5
        K1=[K1 X1(:,i)];
    else
        K2=[K2 X1(:,i)];

    end
end
for i=1:N*3/10              % skup X2
    if pom(2,i)<0.5
        K1=[K1 X2(:,i)];
    else
        K2=[K2 X2(:,i)];
    end
end

        
figure(22); 
plot(K1(1,:),K1(2,:),'o','color',[0.8500 0.3250 0.0980]);
hold on; 
plot(K2(1,:),K2(2,:),'o','color',[0.4660 0.6740 0.1880]);
title('Pocetna klasterizacija (apriorno predznanje) ');
xlabel('X1');
ylabel('X2');
legend('K1','K2');
hold off

%% Klasterizacija sa 2 klase metodom kvadratne dekompozicije

%broj elemenata u klasterima

N1=length(K1(1,:));
N2=length(K2(1,:));

%srednje vrednosti klastera

M1=mean(K1,2);
M2=mean(K2,2);

%kovarijacione matrice klastera

S1=cov(K1');
S2=cov(K2');

%verovatnoce pojavljivanja klasa

P1=length(K1(1,:))/N;
P2=length(K2(1,:))/N;

itermax=100;                           % maksimalan broj iteracija
l=1;
reclust=1;                             % 1 je ako ima reklasterizacije, 0 ako je nema

while (l<itermax) && reclust
    
    K1pom=[];
    K2pom=[];
    
    reclust=0;
    
    for i=1:N1                       % racunamo statisticko rastojanje i trazimo minimum za 1. klasu
        d1=-log(P1)+0.5*log(det(S1))+0.5*(K1(:,i)-M1)'*(S1)^(-1)*(K1(:,i)-M1);
        d2=-log(P2)+0.5*log(det(S2))+0.5*(K1(:,i)-M2)'*(S2)^(-1)*(K1(:,i)-M2);
        dmin=min([d1 d2]);   
        if d1==dmin
            K1pom=[K1pom K1(:,i)];
        else
            K2pom=[K2pom K1(:,i)];
            reclust=1;
        end
    end
    
    for i=1:N2                      % racunamo statisticko rastojanje i trazimo minimum za 2. klasu
        d1=-log(P1)+0.5*log(det(S1))+0.5*(K2(:,i)-M1)'*(S1)^(-1)*(K2(:,i)-M1);
        d2=-log(P2)+0.5*log(det(S2))+0.5*(K2(:,i)-M2)'*(S2)^(-1)*(K2(:,i)-M2);
        dmin=min([d1 d2]);       
        if d1==dmin
            K1pom=[K1pom K2(:,i)];
            reclust=1;
        else
            K2pom=[K2pom K2(:,i)];
        end
    end
    
    clear K1 K2;
    K1=K1pom;
    K2=K2pom;

    figure(23); 
    plot(K1(1,:),K1(2,:),'o','color',[0.8500 0.3250 0.0980]);
    hold on; 
    plot(K2(1,:),K2(2,:),'o','color',[0.4660 0.6740 0.1880]);
    title(['Iteracija broj: ' num2str(l)]);
    pause 
    
    % priprema za narednu iteraciju
    
    N1=length(K1(1,:)); 
    N2=length(K2(1,:));
    
    M1=mean(K1,2);
    M2=mean(K2,2);
    
    S1=cov(K1');
    S2=cov(K2');

    P1=length(K1(1,:))/N;
    P2=length(K2(1,:))/N;
    
    l=l+1;   
end


%% Nasumicna pocetna klasterizacija sa 3 klase

pom=rand(2,N);

K1=[];
K2=[];
K3=[];

for i=1:N                           % skup X1
    if pom(1,i)<1/3
        K1=[K1 X1(:,i)];
    elseif pom(1,i)<2/3
        K2=[K2 X1(:,i)];     
    else 
        K3=[K3 X1(:,i)];
    end
end
for i=1:N                           % skup X2
    if pom(2,i)<1/3
        K1=[K1 X2(:,i)];
    elseif pom(2,i)<2/3
        K2=[K2 X2(:,i)];
    else 
        K3=[K3 X2(:,i)];
    end
end


figure(24); 
plot(K1(1,:),K1(2,:),'o','color',[0.8500 0.3250 0.0980]);
hold on; 
plot(K2(1,:),K2(2,:),'o','color',[0.4660 0.6740 0.1880]); 
plot(K3(1,:),K3(2,:),'o','color',[0 0.4470 0.7410]);
title('Pocetna klasterizacija');
xlabel('X1');
ylabel('X2');
legend('K1','K2','K3');
hold off

%% Klasterizacija sa 3 klase metodom kvadratne dekompozicije

%broj elemenata u klasterima

N1=length(K1(1,:));
N2=length(K2(1,:));
N3=length(K3(1,:));

%srednje vrednosti klastera

M1=mean(K1,2);
M2=mean(K2,2);
M3=mean(K3,2);

%kovarijacione matrice klastera

S1=cov(K1');
S2=cov(K2');
S3=cov(K3');

%verovatnoce pojavljivanja klasa

P1=length(K1(1,:))/N;
P2=length(K2(1,:))/N;
P3=length(K3(1,:))/N;

itermax=100;                        % maksimalan broj iteracija
l=1;
reclust=1;                          % 1 je ako ima reklasterizacije, 0 ako je nema

while (l<itermax) && reclust
    
    K1pom=[];
    K2pom=[];
    K3pom=[];
    
    reclust=0;
    
    for i=1:N1                      % racunamo statisticko rastojanje i trazimo minimum za 1. klasu
        d1=-log(P1)+0.5*log(det(S1))+0.5*(K1(:,i)-M1)'*(S1)^(-1)*(K1(:,i)-M1);
        d2=-log(P2)+0.5*log(det(S2))+0.5*(K1(:,i)-M2)'*(S2)^(-1)*(K1(:,i)-M2);
        d3=-log(P3)+0.5*log(det(S3))+0.5*(K1(:,i)-M3)'*(S3)^(-1)*(K1(:,i)-M3);
        dmin=min([d1 d2 d3]);      
        if d1==dmin
            K1pom=[K1pom K1(:,i)];
        elseif d2==dmin
            K2pom=[K2pom K1(:,i)];
            reclust=1;
        else 
            K3pom=[K3pom K1(:,i)];
            reclust=1;
        end
    end
    
    for i=1:N2                      % racunamo statisticko rastojanje i trazimo minimum za 2. klasu
        d1=-log(P1)+0.5*log(det(S1))+0.5*(K2(:,i)-M1)'*(S1)^(-1)*(K2(:,i)-M1);
        d2=-log(P2)+0.5*log(det(S2))+0.5*(K2(:,i)-M2)'*(S2)^(-1)*(K2(:,i)-M2);
        d3=-log(P3)+0.5*log(det(S3))+0.5*(K2(:,i)-M3)'*(S3)^(-1)*(K2(:,i)-M3);
        dmin=min([d1 d2]);      
        if d1==dmin
            K1pom=[K1pom K2(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K2(:,i)];
        else
            K3pom=[K3pom K2(:,i)];
            reclust=1;
        end
    end
    
    for i=1:N3                      % racunamo statisticko rastojanje i trazimo minimum za 3. klasu
        d1=-log(P1)+0.5*log(det(S1))+0.5*(K3(:,i)-M1)'*(S1)^(-1)*(K3(:,i)-M1);
        d2=-log(P2)+0.5*log(det(S2))+0.5*(K3(:,i)-M2)'*(S2)^(-1)*(K3(:,i)-M2);
        d3=-log(P3)+0.5*log(det(S3))+0.5*(K3(:,i)-M3)'*(S3)^(-1)*(K3(:,i)-M3);
        dmin=min([d1 d2]);       
        if d1==dmin
            K1pom=[K1pom K3(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K3(:,i)];
            reclust=1;
        else
            K3pom=[K3pom K3(:,i)];
        end
    end
    
    clear K1 K2 K3;
    K1=K1pom;
    K2=K2pom;
    K3=K3pom;

    figure(25); 
    plot(K1(1,:),K1(2,:),'o','color',[0.8500 0.3250 0.0980]);
    hold on; 
    plot(K2(1,:),K2(2,:),'o','color',[0.4660 0.6740 0.1880]); 
    plot(K3(1,:),K3(2,:),'o','color',[0 0.4470 0.7410]);
    hold off;
    title(['Iteracija broj: ' num2str(l)]);
    pause 
    
    % priprema za narednu iteraciju
    
    N1=length(K1(1,:)); 
    N2=length(K2(1,:));
    N3=length(K3(1,:));
    
    M1=mean(K1,2);
    M2=mean(K2,2);
    M3=mean(K3,2);
    
    S1=cov(K1');
    S2=cov(K2');
    S3=cov(K3');

    P1=length(K1(1,:))/N;
    P2=length(K2(1,:))/N;
    P3=length(K3(1,:))/N;
    
    l=l+1;   
end





