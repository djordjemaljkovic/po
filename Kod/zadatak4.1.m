%% Domaci zadatak iz prepoznavanja oblika - zadatak 4, zahtev 1

clear;
close all;
clc;

%% Generisanje dvodimenzionalna odbiraka od 4 linearno separabilne klase
N=500;

%Klasa 1

M1 = [-5; 4];
S1 = [3.5 1.1;1.1 2];
[F1,L1]=eig(S1);
T1=F1*L1^(1/2);
X1 = (T1*randn(2,N)+M1*ones(1,N));

%Klasa 2

M2 = [7; 4]; 
S2 = [3 0.8;0.8 0.5];
[F2,L2]=eig(S2);
T2=F2*L2^(1/2);
X2 = (T2*randn(2,N)+M2*ones(1,N));

%Klasa 3

M3 = [-8; -5];
S3 = [2 0.3;0.3 0.5];
[F3,L3]=eig(S3);
T3=F3*L3^(1/2);
X3 = (T3*randn(2,N)+M3*ones(1,N));

%Klasa 4

M4 = [7; -4];
S4 = [4 -0.6;-0.6 0.23];
[F4,L4]=eig(S4);
T4=F4*L4^(1/2);
X4 = (T4*randn(2,N)+M4*ones(1,N));

figure(12);
plot(X1(1,:),X1(2,:),'ro'); 
hold on
plot(X2(1,:),X2(2,:),'go');
plot(X3(1,:),X3(2,:),'bo');
plot(X4(1,:),X4(2,:),'o','color',[0.9290 0.6940 0.1250]);

xlabel('X1');
ylabel('X2');
legend('K1','K2','K3','K4');
title('Odbirci klasa');
hold off

%% C-mean klasterizacija (nasumicna pocetna klasterizacija)

test=rand(4,N);

K1=[];
K2=[];
K3=[];
K4=[];

for i=1:N                       % iteriranje kroz X1
    if test(1,i)<0.25
        K1=[K1 X1(:,i)];
    elseif test(1,i)<0.5
        K2=[K2 X1(:,i)];
    elseif test(1,i)<0.75
        K3=[K3 X1(:,i)];
    else 
        K4=[K4 X1(:,i)];
    end
end

for i=1:N                       % iteriranje kroz X2
    if test(2,i)<0.25
        K1=[K1 X2(:,i)];
    elseif test(2,i)<0.5
        K2=[K2 X2(:,i)];
    elseif test(2,i)<0.75
        K3=[K3 X2(:,i)];
    else 
        K4=[K4 X2(:,i)];
    end
end

for i=1:N                       % iteriranje kroz X3
    if test(3,i)<0.25
        K1=[K1 X3(:,i)];
    elseif test(3,i)<0.5
        K2=[K2 X3(:,i)];
    elseif test(3,i)<0.75
        K3=[K3 X3(:,i)];
    else 
        K4=[K4 X3(:,i)];
    end
end

for i=1:N                        % iteriranje kroz X4
    if test(4,i)<0.25
        K1=[K1 X4(:,i)];
    elseif test(4,i)<0.5
        K2=[K2 X4(:,i)];
    elseif test(4,i)<0.75
        K3=[K3 X4(:,i)];
    else 
        K4=[K4 X4(:,i)];
    end
end
        
figure(13); 
plot(K1(1,:),K1(2,:),'ro');
hold on; 
plot(K2(1,:),K2(2,:),'go'); 
plot(K3(1,:),K3(2,:),'bo');
plot(K4(1,:),K4(2,:),'o','color',[0.9290 0.6940 0.1250]); 
title('Pocetna klasterizacija');
xlabel('X1');
ylabel('X2');
legend('K1','K2','K3','K4');
hold off


%% C-mean klasterizacija (pocetna klasterizacija sa apriornim predznanjem)

test=rand(4,N*3/10);

K1=X1(:,N*3/10+1:end);
K2=X2(:,N*3/10+1:end);
K3=X3(:,N*3/10+1:end);
K4=X4(:,N*3/10+1:end);

for i=1:N*3/10                  % iteriranje kroz X1
    if test(1,i)<0.25
        K1=[K1 X1(:,i)];
    elseif test(1,i)<0.5
        K2=[K2 X1(:,i)];
    elseif test(1,i)<0.75
        K3=[K3 X1(:,i)];
    else 
        K4=[K4 X1(:,i)];
    end
end
for i=1:N*3/10                  % iteriranje kroz X2
    if test(2,i)<0.25
        K1=[K1 X2(:,i)];
    elseif test(2,i)<0.5
        K2=[K2 X2(:,i)];
    elseif test(2,i)<0.75
        K3=[K3 X2(:,i)];
    else 
        K4=[K4 X2(:,i)];
    end
end

for i=1:N*3/10                  % iteriranje kroz X3
    if test(3,i)<0.25
        K1=[K1 X3(:,i)];
    elseif test(3,i)<0.5
        K2=[K2 X3(:,i)];
    elseif test(3,i)<0.75
        K3=[K3 X3(:,i)];
    else 
        K4=[K4 X3(:,i)];
    end
end

for i=1:N*3/10                  % iteriranje kroz X4
    if test(4,i)<0.25
        K1=[K1 X4(:,i)];
    elseif test(4,i)<0.5
        K2=[K2 X4(:,i)];
    elseif test(4,i)<0.75
        K3=[K3 X4(:,i)];
    else 
        K4=[K4 X4(:,i)];
    end
end
        
figure(14); 
plot(K1(1,:),K1(2,:),'ro');
hold on; 
plot(K2(1,:),K2(2,:),'go'); 
plot(K3(1,:),K3(2,:),'bo'); 
plot(K4(1,:),K4(2,:),'o','color',[0.9290 0.6940 0.1250]); 
title('Pocetna klasterizacija sa apriornim predznanjem');
xlabel('X1');
ylabel('X2');
legend('K1','K2','K3','K4');
hold off

%% Klasterizacija 4 klase

%broj elemenata u klasterima

N1=length(K1(1,:));
N2=length(K2(1,:));
N3=length(K3(1,:));
N4=length(K4(1,:));

%srednje vrednosti klastera

M1=mean(K1,2);
M2=mean(K2,2);
M3=mean(K3,2);
M4=mean(K4,2);

itermax=100;               %maks broj iteracija
l=1;
reclust=1;                  % 1 je ako ima reklasterizacije, 0 ako je nema

while (l<itermax) && reclust
    K1pom=[];
    K2pom=[];
    K3pom=[];
    K4pom=[];
    
    reclust=0;
    
    for i=1:N1                  % racunamo euklidsko rastojanje i trazimo minimum za 1. klasu
        d1=sum((K1(:,i)-M1).^2);
        d2=sum((K1(:,i)-M2).^2);
        d3=sum((K1(:,i)-M3).^2);
        d4=sum((K1(:,i)-M4).^2);
        dmin=min([d1 d2 d3 d4]);
        
        if d1==dmin
            K1pom=[K1pom K1(:,i)];
        elseif d2==dmin
            K2pom=[K2pom K1(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K1(:,i)];
            reclust=1;  
        else
            K4pom=[K4pom K1(:,i)];
            reclust=1;  
        end
    end
    
    for i=1:N2                  % racunamo euklidsko rastojanje i trazimo minimum za 2. klasu
        d1=sum((K2(:,i)-M1).^2);
        d2=sum((K2(:,i)-M2).^2);
        d3=sum((K2(:,i)-M3).^2);
        d4=sum((K2(:,i)-M4).^2);
        dmin=min([d1 d2 d3 d4]);
        
        if d1==dmin
            K1pom=[K1pom K2(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K2(:,i)];
        elseif d3==dmin
            K3pom=[K3pom K2(:,i)];
            reclust=1;  
        else
            K4pom=[K4pom K2(:,i)];
            reclust=1;  
        end
    end
    
        for i=1:N3              % racunamo euklidsko rastojanje i trazimo minimum za 3. klasu
        d1=sum((K3(:,i)-M1).^2);
        d2=sum((K3(:,i)-M2).^2);
        d3=sum((K3(:,i)-M3).^2);
        d4=sum((K3(:,i)-M4).^2);
        dmin=min([d1 d2 d3 d4]);
        
        if d1==dmin
            K1pom=[K1pom K3(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K3(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K3(:,i)]; 
        else
            K4pom=[K4pom K3(:,i)];
            reclust=1;  
        end
    end

    for i=1:N4                  % racunamo euklidsko rastojanje i trazimo minimum za 4. klasu
        d1=sum((K4(:,i)-M1).^2);
        d2=sum((K4(:,i)-M2).^2);
        d3=sum((K4(:,i)-M3).^2);
        d4=sum((K4(:,i)-M4).^2);
        dmin=min([d1 d2 d3 d4]);
        
        if d1==dmin
            K1pom=[K1pom K4(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K4(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K4(:,i)];
            reclust=1;
        else
            K4pom=[K4pom K4(:,i)];
        end
    end
    
    clear K1 K2 K3 K4;
    K1=K1pom;
    K2=K2pom;
    K3=K3pom;
    K4=K4pom;
    
    figure(15); 
    plot(K1(1,:),K1(2,:),'ro'); 
    hold on; 
    plot(K2(1,:),K2(2,:),'go');
    plot(K3(1,:),K3(2,:),'bo');
    plot(K4(1,:),K4(2,:),'o','color',[0.9290 0.6940 0.1250]);
    hold off;
    title(['Iteracija broj: ' num2str(l)]);
    pause
    
    % priprema za narednu iteraciju
    
    N1=length(K1(1,:)); 
    N2=length(K2(1,:));
    N3=length(K3(1,:));
    N4=length(K4(1,:));
    
    M1=mean(K1,2);
    M2=mean(K2,2);
    M3=mean(K3,2);
    M4=mean(K4,2);
    
    l=l+1;
  
    
end

%% Klasterizacija kada apriorno ne poznajemo broj klasa ( pretpostavka da postoje 3 klase)


%% C-mean klasterizacija (nasumicna pocetna klasterizacija sa 3 klase)

test=rand(4,N);

K1=[];
K2=[];
K3=[];

for i=1:N                   % skup X1
    if test(1,i)<1/3
        K1=[K1 X1(:,i)];
    elseif test(1,i)<2/3
        K2=[K2 X1(:,i)];     
    else 
        K3=[K3 X1(:,i)];
    end
end
for i=1:N                   % skup X2
    if test(2,i)<1/3
        K1=[K1 X2(:,i)];
    elseif test(2,i)<2/3
        K2=[K2 X2(:,i)];
    else 
        K3=[K3 X2(:,i)];
    end
end

for i=1:N                    % skup X3
    if test(3,i)<1/3
        K1=[K1 X3(:,i)];
    elseif test(3,i)<2/3
        K2=[K2 X3(:,i)];
    else 
        K3=[K3 X3(:,i)];
    end
end

for i=1:N                     % skup X4
    if test(4,i)<1/3
        K1=[K1 X4(:,i)];
    elseif test(4,i)<2/3
        K2=[K2 X4(:,i)];
    else
        K3=[K3 X4(:,i)];
    end
end
        
figure(16); 
plot(K1(1,:),K1(2,:),'ro'); 
hold on; 
plot(K2(1,:),K2(2,:),'go');
plot(K3(1,:),K3(2,:),'bo');
hold off;
title('Pocetna klasterizacija');
xlabel('X1');
ylabel('X2');
legend('K1','K2','K3');
hold off

%% Klasterizacija sa 3 klase

%broj elemenata u klasterima

N1=length(K1(1,:));
N2=length(K2(1,:));
N3=length(K3(1,:));

%srednje vrednosti klastera

M1=mean(K1,2);
M2=mean(K2,2);
M3=mean(K3,2);

itermax=100;                %maksimalni broj iteracija
l=1;
reclust=1;                  % 1 je ako ima reklasterizacije, 0 ako je nema

while (l<itermax) && reclust
    
    K1pom=[];
    K2pom=[];
    K3pom=[];
    
    reclust=0;
    
    for i=1:N1                  % racunamo euklidsko rastojanje i trazimo minimum za 1. klasu
        d1=sum((K1(:,i)-M1).^2);
        d2=sum((K1(:,i)-M2).^2);
        d3=sum((K1(:,i)-M3).^2);
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
    
    for i=1:N2                  % racunamo euklidsko rastojanje i trazimo minimum za 2. klasu
        d1=sum((K2(:,i)-M1).^2);
        d2=sum((K2(:,i)-M2).^2);
        d3=sum((K2(:,i)-M3).^2);
        dmin=min([d1 d2 d3]);
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
    
    for i=1:N3                  % racunamo euklidsko rastojanje i trazimo minimum za 3. klasu
        d1=sum((K3(:,i)-M1).^2);
        d2=sum((K3(:,i)-M2).^2);
        d3=sum((K3(:,i)-M3).^2);
        dmin=min([d1 d2 d3]);   
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
    
    figure(17); 
    plot(K1(1,:),K1(2,:),'ro'); 
    hold on; 
    plot(K2(1,:),K2(2,:),'go');
    plot(K3(1,:),K3(2,:),'bo');
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

    l=l+1;
end

%% Klasterizacija kada apriorno ne poznajemo broj klasa ( pretpostavka da postoji 5 klasa) 


%% C-mean klasterizacija (nasumicna pocetna klasterizacija sa 5 klasa)

test=rand(4,N);

K1=[];
K2=[];
K3=[];
K4=[];
K5=[];

for i=1:N                   % skup X1
    if test(1,i)<1/5
        K1=[K1 X1(:,i)];
    elseif test(1,i)<2/5
        K2=[K2 X1(:,i)];
    elseif test(1,i)<3/5
        K3=[K3 X1(:,i)];
    elseif  test(1,i)<4/5
        K4=[K4 X1(:,i)];
    else
        K5=[K5 X1(:,i)];
    end
end
for i=1:N                   % skup X2
    if test(2,i)<1/5
        K1=[K1 X2(:,i)];
    elseif test(2,i)<2/5
        K2=[K2 X2(:,i)];
    elseif test(2,i)<3/5
        K3=[K3 X2(:,i)];
    elseif  test(2,i)<4/5
        K4=[K4 X2(:,i)];
    else
        K5=[K5 X2(:,i)];
    end
end

for i=1:N                   % skup X3
    if test(3,i)<1/5
        K1=[K1 X3(:,i)];
    elseif test(3,i)<2/5
        K2=[K2 X3(:,i)];
    elseif test(3,i)<3/5
        K3=[K3 X3(:,i)];
    elseif  test(3,i)<4/5
        K4=[K4 X3(:,i)];
    else
        K5=[K5 X3(:,i)];
    end
end

for i=1:N                   % skup X4
    if test(4,i)<1/5
        K1=[K1 X4(:,i)];
    elseif test(4,i)<2/5
        K2=[K2 X4(:,i)];
    elseif test(4,i)<3/5
        K3=[K3 X4(:,i)];
    elseif  test(4,i)<4/5
        K4=[K4 X4(:,i)];
    else
        K5=[K5 X4(:,i)];
    end
end
        
figure(18); 
plot(K1(1,:),K1(2,:),'ro'); 
hold on; 
plot(K2(1,:),K2(2,:),'go');
plot(K3(1,:),K3(2,:),'bo');
plot(K4(1,:),K4(2,:),'o', 'color',[0.9290 0.6940 0.1250] ); 
plot(K5(1,:),K5(2,:),'mo'); 

title('Pocetna klasterizacija');
xlabel('X1');
ylabel('X2');
legend('K1','K2','K3','K4','K5');
hold off

%% Klasterizacija sa 5 klasa

%broj elemenata u klasterima

N1=length(K1(1,:));
N2=length(K2(1,:));
N3=length(K3(1,:));
N4=length(K4(1,:));
N5=length(K5(1,:));

%srednje vrednosti klastera

M1=mean(K1,2);
M2=mean(K2,2);
M3=mean(K3,2);
M4=mean(K4,2);
M5=mean(K5,2);

itermax=100;                    %maksimalan broj iteracija
l=1;
reclust=1;                      % 1 je ako ima reklasterizacije, 0 ako je nema

while (l<itermax) && reclust
    
    K1pom=[];
    K2pom=[];
    K3pom=[];
    K4pom=[];
    K5pom=[];
    
    reclust=0;
    
    for i=1:N1                  % racunamo euklidsko rastojanje i trazimo minimum za 1. klasu
        d1=sum((K1(:,i)-M1).^2);
        d2=sum((K1(:,i)-M2).^2);
        d3=sum((K1(:,i)-M3).^2);
        d4=sum((K1(:,i)-M4).^2);
        d5=sum((K1(:,i)-M5).^2);
        dmin=min([d1 d2 d3 d4 d5]);        
        if d1==dmin
            K1pom=[K1pom K1(:,i)];
        elseif d2==dmin
            K2pom=[K2pom K1(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K1(:,i)];
            reclust=1;  
        elseif d4==dmin
            K4pom=[K4pom K1(:,i)];
            reclust=1;  
        else
            K5pom=[K5pom K1(:,i)];
            reclust=1;
        end
    end
    
    for i=1:N2                  % racunamo euklidsko rastojanje i trazimo minimum za 2. klasu
        d1=sum((K2(:,i)-M1).^2);
        d2=sum((K2(:,i)-M2).^2);
        d3=sum((K2(:,i)-M3).^2);
        d4=sum((K2(:,i)-M4).^2);
        d5=sum((K2(:,i)-M5).^2);
        dmin=min([d1 d2 d3 d4 d5]);       
        if d1==dmin
            K1pom=[K1pom K2(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K2(:,i)];
        elseif d3==dmin
            K3pom=[K3pom K2(:,i)];
            reclust=1;  
        elseif d4==dmin
            K4pom=[K4pom K2(:,i)];
            reclust=1;  
        else
            K5pom=[K5pom K2(:,i)];
            reclust=1;
        end
        
    end
    
    for i=1:N3                  % racunamo euklidsko rastojanje i trazimo minimum za 3. klasu
        d1=sum((K3(:,i)-M1).^2);
        d2=sum((K3(:,i)-M2).^2);
        d3=sum((K3(:,i)-M3).^2);
        d4=sum((K3(:,i)-M4).^2);
        d5=sum((K3(:,i)-M5).^2);
        dmin=min([d1 d2 d3 d4 d5]);      
        if d1==dmin
            K1pom=[K1pom K3(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K3(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K3(:,i)]; 
        elseif d4==dmin
            K4pom=[K4pom K3(:,i)];
            reclust=1;
        else
            K5pom=[K5pom K3(:,i)];
            reclust=1;
        end
    end

    for i=1:N4                  % racunamo euklidsko rastojanje i trazimo minimum za 4. klasu
        d1=sum((K4(:,i)-M1).^2);
        d2=sum((K4(:,i)-M2).^2);
        d3=sum((K4(:,i)-M3).^2);
        d4=sum((K4(:,i)-M4).^2);
        d5=sum((K4(:,i)-M5).^2);
        dmin=min([d1 d2 d3 d4 d5]);     
        if d1==dmin
            K1pom=[K1pom K4(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K4(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K4(:,i)];
            reclust=1;
        elseif d4==dmin
            K4pom=[K4pom K4(:,i)];
        else
            K5pom=[K5pom K4(:,i)];
            reclust=1;
        end
    end
    
    for i=1:N5                  % racunamo euklidsko rastojanje i trazimo minimum za 5. klasu
        d1=sum((K5(:,i)-M1).^2);
        d2=sum((K5(:,i)-M2).^2);
        d3=sum((K5(:,i)-M3).^2);
        d4=sum((K5(:,i)-M4).^2);
        d5=sum((K5(:,i)-M5).^2);
        dmin=min([d1 d2 d3 d4 d5]);        
        if d1==dmin
            K1pom=[K1pom K5(:,i)];
            reclust=1;
        elseif d2==dmin
            K2pom=[K2pom K5(:,i)];
            reclust=1;
        elseif d3==dmin
            K3pom=[K3pom K5(:,i)];
            reclust=1;
        elseif d4==dmin
            K4pom=[K4pom K5(:,i)];
            reclust=1;
        else
            K5pom=[K5pom K5(:,i)];
        end
    end
    
    
    clear K1 K2 K3 K4 K5;
    K1=K1pom;
    K2=K2pom;
    K3=K3pom;
    K4=K4pom;
    K5=K5pom;
    
    figure(19); 
    plot(K1(1,:),K1(2,:),'ro'); 
    hold on; 
    plot(K2(1,:),K2(2,:),'go');
    plot(K3(1,:),K3(2,:),'bo');
    plot(K4(1,:),K4(2,:),'o', 'color',[0.9290 0.6940 0.1250] ); 
    plot(K5(1,:),K5(2,:),'mo'); 
    hold off
    title(['Iteracija broj: ' num2str(l)]);
    pause 
    
    % priprema za narednu iteraciju
    
    N1=length(K1(1,:)); 
    N2=length(K2(1,:));
    N3=length(K3(1,:));
    N4=length(K4(1,:));
    N5=length(K5(1,:));
    
    M1=mean(K1,2);
    M2=mean(K2,2);
    M3=mean(K3,2);
    M4=mean(K4,2);
    M5=mean(K5,2);
    
    l=l+1;
end
