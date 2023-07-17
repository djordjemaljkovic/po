
%% Doma?i zadatak iz prepoznavanja oblika - zadatak 1

clear
close all
clc
cd BazaSlova;
%% ucitavanje svih slova A


x=dir('bazaA*.bmp'); 
for i=1:max(size(x))
    X=imread(x(i).name);
    P1(:,i)=Obelezja(X);    %obelezja za A
end

%% ucitavanje svih slova E


x=dir('bazaE*.bmp');
for i=1:max(size(x))
    X=imread(x(i).name); 
    P2(:,i)=Obelezja(X); %obelezja za E
end



%% ucitavanje svih slova I


x=dir('bazaI*.bmp');
for i=1:max(size(x))
    X=imread(x(i).name); 
    P3(:,i)=Obelezja(X);
end


%% ucitavanje svih slova O

x=dir('bazaO*.bmp');
for i=1:max(size(x))
    X=imread(x(i).name); 
    P4(:,i)=Obelezja(X);
end


%% ucitavanje svih slova U

x=dir('bazaU*.bmp');
for i=1:max(size(x))
    X=imread(x(i).name); 
    P5(:,i)=Obelezja(X);
end


%% Iscrtavanje obelezja

figure(4)
plot3(P1(1,:),P1(2,:),P1(3,:),'ro');
hold on;
plot3(P2(1,:),P2(2,:),P2(3,:),'bx');
plot3(P3(1,:),P3(2,:),P3(3,:),'g*');
plot3(P4(1,:),P4(2,:),P4(3,:),'mx');
plot3(P5(1,:),P5(2,:),P5(3,:),'yv');
hold off;
legend('A','E','I','O','U');
xlabel('1. obelezje');
ylabel('2. obelezje');
zlabel('3. obelezje');

%% Klasifikacija


N=120;
No=100;                             % obucavajuci skup i testirajuci skup
O1=P1(:,1:No);T1=P1(:,No+1:N);
O2=P2(:,1:No);T2=P2(:,No+1:N);
O3=P3(:,1:No);T3=P3(:,No+1:N);
O4=P4(:,1:No);T4=P4(:,No+1:N);
O5=P5(:,1:No);T5=P5(:,No+1:N);

M1=mean(O1,2); S1=cov(O1');         %P1=P2=P3=P4=P5, raspodele su normalne
M2=mean(O2,2); S2=cov(O2');
M3=mean(O3,2); S3=cov(O3');
M4=mean(O4,2); S4=cov(O4');
M5=mean(O5,2); S5=cov(O5');

% Testiranje klasifikatora: za svaku klasu, testiramo svaki odbirak i belezimo da li je odluka klasifikatora ispravna ili ne

Mk=zeros(5);                    % konfuziona matrica
for klasa=1:5
    disp(num2str(klasa));
    if klasa==1
        T=T1;
    elseif klasa==2
        T=T2;
    elseif klasa==3
        T=T3;
    elseif klasa==4
        T=T4;
    else
        T=T5;
    end
    for i=1:N-No           
        p=T(:,i);
        f1=1/2/pi/det(S1)^0.5*exp(-0.5*(p-M1)'*inv(S1)*(p-M1));
        f2=1/2/pi/det(S2)^0.5*exp(-0.5*(p-M2)'*inv(S2)*(p-M2));
        f3=1/2/pi/det(S3)^0.5*exp(-0.5*(p-M3)'*inv(S3)*(p-M3));
        f4=1/2/pi/det(S4)^0.5*exp(-0.5*(p-M4)'*inv(S4)*(p-M4));
        f5=1/2/pi/det(S5)^0.5*exp(-0.5*(p-M5)'*inv(S5)*(p-M5));
        %odluku donosi maksimalni element
        m=max([f1 f2 f3 f4 f5]); 	
        if m==f1                            % klasifikator doneo odluku da je klasa 1
            Mk(klasa,1)=Mk(klasa,1)+1;
            disp('Odluka je 1');
        elseif m==f2                        % klasifikator doneo odluku da je klasa 2
            Mk(klasa,2)=Mk(klasa,2)+1;
            disp('Odluka je 2');
        elseif m==f3                        % klasifikator doneo odluku da je klasa 3
            Mk(klasa,3)=Mk(klasa,3)+1;
            disp('Odluka je 3');
        elseif m==f4                        % klasifikator doneo odluku da je klasa 4
            Mk(klasa,4)=Mk(klasa,4)+1;
            disp('Odluka je 4');
        else                                % klasifikator doneo odluku da je klasa 5
            Mk(klasa,5)=Mk(klasa,5)+1;
            disp('Odluka je 5');
        end
    end
end

Mk
disp(['Greska iznosi: ' num2str((sum(sum(Mk))-trace(Mk))/sum(sum(Mk)))]);