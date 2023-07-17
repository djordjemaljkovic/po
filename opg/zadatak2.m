%% Obrada i prepoznavanje govora zadatak 2

clear
close all; 
clc;

%% Defining constants
fs = 8000;         %ucestanost odabiranja
T = 1/fs;          %period odabiranja
duration = 12;     %duzina govorne sekvence
N = duration*fs;   %broj uzoraka
nbits=16;          %broj bitova za kodiranje
nchans = 1;        %broj govornih kanala

%% Snimanje govornog signala

 x = audiorecorder(fs, nbits, nchans);
 disp('Pocetak')
 recordblocking(x, duration);
 disp('Kraj')

%% Cuvanje audio sekvence

rec_name = 'zedatak3.wav';
y = getaudiodata(x);
audiowrite(rec_name, y, fs);
[x, Fs] = audioread(rec_name);
N = length(x);
t = 0:1/Fs:(N-1)/fs;
sound(x, fs);

%% mi companding


Xmax = max(abs(x));
fig=1;

for mi=[100 500]
    xquant_fcn=Xmax*log10(1+mi*abs(x)/Xmax)/(log10(1+mi)).*sign(x); 
    j=1;
    for b=[4 8 12]
        M = 2^b;             %kvantizacioni nivoi
        d = 2*Xmax/M;        %razlika susednih nivoa kvantizatora
        
        xquant_mi=round(xquant_fcn/d)*d; 
        xquant_mi(xquant_fcn>(M-1)*d/2)=(M/2-1)*d;
        xq_mid = 1/mi*sign(xquant_mi).*((1+mi).^(abs(xquant_mi)/Xmax)-1)*Xmax; 
        
        %sound(xq_mid,fs);
        figure(fig)
        subplot(3,1,j)
        plot(x,xquant_mi,'*');   
        ylabel(['b=' + string(b)]);
        
        figure(fig+1)   
        subplot(3,1,j)
        plot(t,x,'b');
        hold on;
        plot(t, xq_mid,'r--');
        ylabel(['b=' + string(b)]);
        j = j+1;
    end
    fig=fig+2;
end

figure(1);
subplot(3,1,1);
title('Karakteristike kvantizatora za \mu=100')
subplot(3,1,3); 
xlabel('t[s]');
figure(2);
subplot(3,1,1);
title('Originalan i kvantiziran signal za \mu=100')
subplot(3,1,3);
xlabel('t[s]');
figure(3);
subplot(3,1,1);
title('Karakteristike kvantizatora za \mu=500')
subplot(3,1,3); 
xlabel('t[s]');
figure(4);
subplot(3,1,1);
title('Originalan i kvantiziran signal za \mu=500')
subplot(3,1,3);
xlabel('t[s]');

%% mi companding - SNR

[x, Fs] = audioread(rec_name);
Xmax = max(abs(x));
fig=5;

for b = [4 8 12]
    M = 2^b;                   %kvantizacioni nivoi
    d = 2*Xmax/M;              %razlika izmedju susednih nivoa
    for mi=[100 500]
        slabljenje = 0.01:0.01:1;
        xvar = zeros(size(slabljenje));
        SNR_mi = zeros(size(slabljenje));
        for i=1:length(slabljenje)
            x1=x*slabljenje(i);
            xvar(i)=var(x1);
            xquant_fcn=Xmax*log10(1+mi*abs(x1)/Xmax)/(log10(1+mi)).*sign(x1); 
            xquant_mi=round(xquant_fcn/d)*d; 
            xquant_mi(xquant_fcn>(M-1)*d/2)=(M/2-1)*d;
            x2 = 1/mi*sign(xquant_mi).*((1+mi).^(abs(xquant_mi)/Xmax)-1)*Xmax; 
            SNR_mi(i) = var(x1)/var(x1-x2);
        end
        figure(fig);
        semilogx(Xmax./sqrt(xvar), 10*log10(SNR_mi),'b');
        hold on
        semilogx(Xmax./sqrt(xvar), 4.77+6*b-20*log10(log(1+mi))-10*log10(1+sqrt(2)*Xmax./(mi.*sqrt(xvar))+(Xmax./(mi.*sqrt(xvar))).^2),'--r');
        legend('Eksperimentalno SNR','Teoretsko SNR');
        ylabel('SNR[dB]');
    end
    hold off
    title(['SNR za kvantizator od ' + string(b) + ' bita'])
    fig=fig+1;
end

%% Delta-kvantizator

fig = 8;
[x, Fs] = audioread(rec_name);
N = length(x);
t = 0:1/Fs:(N-1)/fs;

delta = [0.08 0.001 0.0105];
xmean = mean(x);

for Q = delta
    
    d = zeros(1,N);         % inkrement
    d(1) = x(1);
    c = zeros(1,N);         % kodirana rec
    xx = zeros(1,N);        % rekonstruisan signal
    xx(1) = xmean+Q;

   
    for i = 2:N
        d(i) = x(i) - xx(i-1);   %predikcija
        if d(i) > 0
            c(i) = 0;
            xx(i) = xx(i-1)+Q;
        else
            c(i) = 1;
            xx(i) = xx(i-1)-Q;
        end
    end

    figure(fig)
    histogram(d);
    title(['Korak kvantizacije = ' + string(Q)]);
    if Q == delta(3)
        figure(fig+1)
        plot(t,x,'*',t,xx,'x');
        legend('Originalni signal', 'Rekonstruisan signal');
        title(['Korak kvantizacije = ' + string(Q)]);
    end
    fig = fig+1;
    
end
%% Sound check
sound(xx,Fs);