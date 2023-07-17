function P = Obelezja(X)        %P je vektor gde ce se smestati obelezja
%% Ucitavanje slika
Y=double(X);
 figure(1);
 imshow(X);
 figure(2);
 imshow(Y/255);


%% Binarizacija

T=0.8;
Y(Y<T*max(max(Y)))=0;               %80 procenata slike
Y(Y>=T*max(max(Y)))=255;            %20 procenata slike 

%% Predobrada slike - uklanjanje okvira

%na nekim slikama pre crnog okvira postoji i beli niz pa se okvir zbog toga
%ne otkloni. Potrebno je i to srediti.

    [redovi,kolone]=size(Y);
    
    %otkljanjanje belog okvira
    
    gore=1; 
    while (gore<redovi) && (sum(Y(gore,:))/kolone >252) || ((sum(Y(gore,:))/kolone<252) && (sum(Y(gore+1,:))/kolone>252)) || ((sum(Y(gore,:))/kolone<252) && (sum(Y(gore+15,:))/kolone<120))
        gore=gore+1;
    end

    dole=redovi;
    while (dole>1) && (sum(Y(dole,:))/kolone>252) || ((sum(Y(dole,:))/kolone<252) && (sum(Y(dole-1,:))/kolone>252))|| ((sum(Y(dole,:))/kolone<252) && (sum(Y(dole-5,:))/kolone<120))
        dole=dole-1;
    end

    levo=1;
    while (levo<kolone) && (sum(Y(:,levo))/redovi>252) || ((sum(Y(:,levo))/redovi<252) && (sum(Y(:,levo+1))/redovi>252)) || ((sum(Y(:,levo))/redovi<252) && (sum(Y(levo+5,:))/redovi<120))
        levo=levo+1;
    end

    desno=kolone;
    while (desno>1) && (sum(Y(:,desno))/redovi>252) || ((sum(Y(:,desno))/redovi<252) && (sum(Y(:,desno-1))/redovi>252)) ||  ((sum(Y(:,desno))/redovi<252) && (sum(Y(desno-5,:))/redovi<120))
        desno=desno-1;
    end
    
    X=X(gore:dole,levo:desno);
    Y=Y(gore:dole,levo:desno);
    [redovi,kolone]=size(Y);
    
    % otklanjanje okvira
    
    gore=1;
    while (gore<redovi) && (sum(Y(gore,:))/kolone<200)
        gore=gore+1;
    end

    dole=redovi;
    while (dole>1) && (sum(Y(dole,:))/kolone<200)
        dole=dole-1;
    end

    levo=1;
    while (levo<kolone) && (sum(Y(:,levo))/redovi<200)
        levo=levo+1;
    end

    desno=kolone;
    while (desno>1) && (sum(Y(:,desno))/redovi<200)
        desno=desno-1;
    end
    
    X=X(gore:dole,levo:desno);
    figure(2);
    imshow(X);
    Y=Y(gore:dole,levo:desno);
    [redovi,kolone]=size(X);

%% Predobrada slike - isecanje belih segmenata

    % isecanje belih segmenata
    
    gore=1;
    while (gore<redovi) && (sum(Y(gore,:))/kolone>252) || ((sum(Y(gore,:))/kolone<252) && (sum(Y(gore+1,:))/kolone>252))
        gore=gore+1;
    end

    dole=redovi;
    while (dole>1) && (sum(Y(dole,:))/kolone>252) ||((sum(Y(dole,:))/kolone<252) && (sum(Y(dole-1,:))/kolone>252))
        dole=dole-1;
    end

    levo=1;
    while (levo<kolone) && (sum(Y(:,levo))/redovi>252) || ((sum(Y(:,levo))/redovi<252) && (sum(Y(:,levo+1))/redovi>252))
        levo=levo+1;
    end

    desno=kolone;
    while (desno>1) && (sum(Y(:,desno))/redovi>252) || ((sum(Y(:,desno))/redovi<252) && (sum(Y(:,desno-1))/redovi>252))
        desno=desno-1;
    end

    X=X(gore:dole,levo:desno);
    Y=Y(gore:dole,levo:desno);

    figure(3);
    imshow(Y/255);
    [redovi,kolone]=size(X);


%% obelezja  

P(1,1) =  mean(mean(Y(round(1/3*redovi):round(2/3*redovi),round(1/3*kolone):round(2/3*kolone))));       %srednji kvadrat
P(2,1) = (mean(mean(Y(1:round(1/5*redovi),round(1/3*kolone):round(2/3*kolone))))+mean(mean(Y(round(1/2*redovi):round(3/4*redovi),1:round(1/2*kolone)))))/2;             %gornji srednji kvadrat plus levi srednji sa dva
P(3,1) = (mean(mean(Y(1:round(1/7*redovi),round(1/3*kolone):round(2/3*kolone))))+mean(mean(Y(round(6/7*redovi):redovi,round(1/3*kolone):round(2/3*kolone)))))/2;        %sredina zbira gornjeg i donjeg srednjeg

