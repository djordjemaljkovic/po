%% Obrada i prepoznavanje govora zadatak 1

clear;
close all; 
clc;

%% Definisanje konstanti govornog signala
fs = 8000;         %ucestanost odabiranja
T = 1/fs;          %period odabiranja
duration = 12;     %duzina govorne sekvence
N = duration*fs;   %broj uzoraka
nbits=16;          %broj bitova za kodiranje
nchans = 1;        %broj govornih kanala

%% Snimanje i cuvanje govorne sekvence
 x = audiorecorder(fs, nbits, nchans);
 disp('Pocetak')
 recordblocking(x, duration);
 disp('Kraj')
 rec_name = 'zadatak1.wav';
 y = getaudiodata(x);
 audiowrite(rec_name, y, fs);

%% Ucitavanje govorne sekvence
[y, fs] = audioread(rec_name);
t = 0:1/fs:(length(y)-1)/fs;
sound(y, fs);

%% Iscrtavanje govorne sekvence
 figure(1)
 plot(t, y); 
 title('Prva govorna sekvenca'); 
 xlabel('t[s]'); 
 ylabel('y(t)');

%% Filtriranje signala
nord = 6; 
Wn = [300 3400]/(fs/2); 
[b, a] = butter(nord, Wn, 'bandpass'); 
yf = filter(b, a, y);

figure(2)
plot(t, y); 
hold on;
plot(t, yf); 
title('Audio signal'); 
xlabel('t[s]'); 
ylabel('y(t)'); 
legend('original', 'filtered');
    
figure(3)
plot(abs(fftshift(fft(y)))); 
hold on;
plot(abs(fftshift(fft(yf)))); 
hold off;
title('Signal spectre'); legend('original', 'filtered');
    
%% Kratkovremenska energija i ZCR
wlength = fs*30e-3+1;       %duzina prozora
window = hamming(wlength);     %Hammingov prozor
E = zeros(size(yf));        %kratkovremenska energija
Z = zeros(size(yf));        %zero-crossing rate

for i = (wlength+1)/2: length(yf)-wlength/2
    range = (i-(wlength-1)/2) : (i+wlength/2);
    E(i) = sum(yf(range).^2.*window);                             %kratkovremenska energija
    Z(i) = sum(abs(sign(yf(range+1))-sign(yf(range))).*window);   %zero-crossng rate
end

Z = Z/(2*wlength);

figure(4);
grafik = plotyy(t,yf,t,E);
title('Talasni oblik i kratkovremenska energija');
xlabel('t[s]');
ylabel(grafik(1), 'Govorni signal');
ylabel(grafik(2), 'Energija');

figure(5); 
grafik = plotyy(t,yf,t,Z);
title('Talasni oblik i ZCR');
xlabel('t[s]');
ylabel(grafik(1), 'Govorni signal');
ylabel(grafik(2), 'ZCR');

%% Segmentacija govornog signala
GEP = 0.08*max(E);
DEP = 0.0014*max(E);

word_beginning = [];
word_end = [];


for i=2:length(E)
    if (E(i-1)<GEP)&&(E(i)>=GEP)
        word_beginning = [word_beginning i];
    end
    if (E(i-1)>GEP)&&(E(i)<=GEP)
        word_end = [word_end i];
    end
end

beginning = word_beginning;
ending = word_end;


rec = zeros(length(E),1);
for i=1:length(beginning)
   rec(beginning(i):ending(i))=max(E)*ones(ending(i)-beginning(i)+1,1); 
end

figure(7)
plot(t,E,'b',t,rec,'r');
    

for i = 1:length(beginning)
   temp = beginning(i);
   while(E(temp)>DEP)
       temp = temp-1;
   end
   beginning(i) = temp;
end

for i = 1:length(ending)
   temp = ending(i);
   while(E(temp)>DEP)
       temp = temp+1;
   end
   ending(i) = temp;
end


clear beginning1 end1
beginning1(1) = beginning(1);
k=1;
for i=2:length(beginning)
   if beginning(i)~= beginning1(k)
       k=k+1;
       beginning1(k) = beginning(i);
   end
end

ending1(1) = ending(1);
k=1;
for i=2:length(ending)
   if ending(i)~= ending1(k)
       k=k+1;
       ending1(k) = ending(i);
   end
end
    
clear rec beginning end
beginning = beginning1;
ending = ending1;
    
rec = zeros(length(yf),1);
for i=1:length(beginning)
   rec(beginning(i):ending(i))=max(E)*ones(ending(i)-beginning(i)+1,1); 
end

hold on;
plot(t,rec/2,'g');
hold off;
legend('Govorni signal', 'Stara procena reci', 'Nova procena reci');
title('Segmentacija govornog signala'); 
xlabel('t[s]');

%% ZCR

z1 = Z(beginning(1)+100:ending(1)-100);    % rec
z2 = Z(2.8*fs:3.2*fs);                     % tisina
h = imhist([z1;z2],40);
Threshold = otsuthresh(h);

figure(8);
hold all;
histogram(z1,20);
histogram(z2,20);
plot([Threshold,Threshold],[0,1000],'r','LineWidth',2);
title('ZCR histogram');
%% Segmentacija govornog signala pomocu ZCR


for i = 1:length(beginning)
    temp = beginning(i);
    count = 0;
    over = false;
    for j = 1:25
        k = beginning(i) - j;
        if k<1
            break;
        end
        
        if over == false && Z(k)<=Threshold
            temp = k;
        elseif over == false && Z(k)>Threshold
            over = true;
        end
        
        if sign(Threshold - Z(k)) ~= sign(Threshold - Z(k+1))
            count = count + 1;
        end     
    end
    if count>=3
        beginning(i) = temp;
    end
end

for i = 1:length(ending)
    temp = ending(i);
    count = 0;
    over = false;
    for j = 1:25
        k = ending(i) + j;
        if k>length(Z)
            break;
        end
        
        if over == false && Z(k)<=Threshold
            temp = k;
        elseif over == false && Z(k)>Threshold
            over = true;
        end
        
        if sign(Threshold - Z(k)) ~= sign(Threshold - Z(k-1))
            count = count + 1;
        end     
    end
   
    if count>=3
        beginning(i) = temp;
    end
end

clear beginning1 ennd1
beginning1(1) = beginning(1);
k=1;
for i=2:length(beginning)
   if beginning(i)~= beginning1(k)
       k=k+1;
       beginning1(k) = beginning(i);
   end
end

ending1(1) = ending(1);
k=1;
for i=2:length(ending)
   if ending(i)~= ending1(k)
       k=k+1;
       ending1(k) = ending(i);
   end
end
    
clear rec beginning end
beginning = beginning1;
ending = ending1;
    
rec = zeros(length(yf),1);
for i=1:length(beginning)
   rec(beginning(i):ending(i))=max(E)*ones(ending(i)-beginning(i)+1,1); 
end

figure(9);
hold on;
plot(t,E,t,rec/2,'g');
hold off;
legend('Govorni signal', 'Konacna segmentacija');
title('Segmentisane reci'); 
xlabel('t[s]');
    
%% Reprodukcija reci

 for i=1:length(beginning)
    sound(y(beginning(i):ending(i)),fs);
    pause
 end

%% Pitch perioda

clear; 
close all; 
clc; 

%% Definisanje konstanti

fs = 8000;         %ucestanost odabiranja
T = 1/fs;          %perioda odabiranja
duration = 10;      %duzina snimanja
N = duration*fs;   %broj odbiraka
nbits=16;          %broj bita za kodiranje
nchans = 1;        %broj kanala

%% Snimanje glasa

x = audiorecorder(fs, nbits, nchans);
disp('Pocetak')
recordblocking(x, duration);
disp('Kraj')
rec_name = 'zadatak2.wav';
y = getaudiodata(x);
audiowrite(rec_name, y, fs);

%% Ucitavanje govornog signala

[y, fs] = audioread(rec_name);
%y = y(0.5*fs:3*fs);
N = length(y);
t = 0:1/fs:(length(y)-1)/fs;
sound(y, fs);
    
%% Filtriranje signala

nord = 6; 
Wn = 300/(fs/2); 
[b, a] = butter(nord, Wn); 
yf = filter(b, a, y);

figure(10)
plot(t, y);
hold on;
plot(t, yf); 
title('Govorni signal'); 
xlabel('t[s]'); 
ylabel('y(t)'); 
legend('originalni', 'filtriran');
    
figure(11)
plot(abs(fftshift(fft(y)))); 
hold on;
plot(abs(fftshift(fft(yf))));
hold off;
title('Spektar govornog signala'); 
legend('originalni', 'filtriran');
    
%% Formiranje impulsa

y = yf;
maxs = zeros(1, N);
mins = zeros(1, N);
maxindex = [];
minindex = [];
 
m1 = zeros(1,N);
m2 = zeros(1,N);
m3 = zeros(1,N);
m4 = zeros(1,N);
m5 = zeros(1,N);
m6 = zeros(1,N);


for i = 2:N-1
    if y(i) > y(i-1) && y(i) > y(i+1)
        maxs(i) = y(i); 
        maxindex = [maxindex i]; 
        m1(i) = max(0, y(i)); 
    end
    if y(i) < y(i-1) && y(i) < y(i+1)
        mins(i) = y(i); 
        minindex = [minindex i]; 
        m4(i) = max(0, -y(i));
    end
end


maxp = 0; 
for i = maxindex
    if isempty(find(minindex < i, 1, 'last'))
        minp = 0; 
    else
        idx = find(minindex < i, 1, 'last'); 
        minp = mins(minindex(idx)); 
    end
    m2(i) = max(0, maxs(i) - minp); 
    m3(i) = max(0, maxs(i) - maxp); 
    maxp = maxs(i); 
end

minp = 0; 
for i = minindex
    if isempty(find(maxindex < i, 1, 'last'))
        maxp = 0; 
    else
        idx = find(maxindex < i, 1, 'last'); 
        maxp = maxs(maxindex(idx)); 
    end
    m5(i) = max(0, -(mins(i) - maxp));
    m6(i) = max(0, -(mins(i) - minp));
    minp = mins(i); 
end

n=1:300;
figure(12)
subplot(3,2,1)
stem(t(n), m1(n)); 
title('m1(t)');
xlabel('t(s)');
subplot(3,2,2) 
stem(t(n), m2(n));
title('m2(t)'); 
xlabel('t(s)');
subplot(3,2,3) 
stem(t(n), m3(n));
title('m3(t)'); 
xlabel('t(s)');
subplot(3,2,4) 
stem(t(n), m4(n)); 
title('m4(t)');
xlabel('t(s)');
subplot(3,2,5) 
stem(t(n), m5(n));
title('m5(t)'); 
xlabel('t(s)');
subplot(3,2,6) 
stem(t(n), m6(n));
title('m6(t)');
xlabel('t(s)');

%% Procena Tp

window = round(fs*20e-3);
NN = floor(N/(window/2));
lambda = 120/fs;            %parametar eksponencijalnog opadanja
tau = round(fs*4e-3);       %"blanking" period  


pt1 = zeros(1, NN);
pt2 = zeros(1, NN);
pt3 = zeros(1, NN); 
pt4 = zeros(1, NN); 
pt5 = zeros(1, NN);
pt6 = zeros(1, NN); 

pt = zeros(1, NN);
k = 2;
for k_win = 1:window/2:N-window+1
    % Prvi estimator
    x = m1(k_win:k_win+window-1);
    pt1(k) = estimate(x, lambda, tau, window, fs);
    % Drugi estimator
    x = m2(k_win:k_win+window-1);
    pt2(k) = estimate(x, lambda, tau, window, fs);
    % Treci estimator
    x = m3(k_win:k_win+window-1);
    pt3(k) = estimate(x, lambda, tau, window, fs);
    % Cetvrti estimator
    x = m4(k_win:k_win+window-1);
    pt4(k) = estimate(x, lambda, tau, window, fs);
    % Peti estimator
    x = m5(k_win:k_win+window-1);
    pt5(k) = estimate(x, lambda, tau, window, fs);
    % Sesti estimator
    x = m6(k_win:k_win+window-1);
    pt6(k) = estimate(x, lambda, tau, window, fs);

    pt(k) = nanmedian([pt1(k) pt2(k) pt3(k) pt4(k) pt5(k) pt6(k) pt(k-1)]);
    k = k + 1; 
end

figure(13)
hold on;
plot(1./pt1); 
plot(1./pt2); 
plot(1./pt3); 
plot(1./pt4); 
plot(1./pt5);
plot(1./pt6); 
plot(1./pt, 'Linewidth', 2); 
title('Procena pitch periode'); 
xlabel('t[s]');
legend('1', '2', '3', '4', '5', '6', 'f_p')
ylim([0,300]);
    
%% Finalna procena Tp

Tp = []; 
for l = 4:length(pt)-3
    Tp = [Tp median(pt(l-3:l+3))]; 
end

figure(14)
plot(1./Tp); 
title('Finalna procena Pitch periode');
xlabel('t[s]');
   
Fp = median(1./Tp);
disp(['Dobijena pitch frekvencija je: ' + string(Fp) + ' Hz']);


%% Odsecanje
start_time = find(y > 0.002, 1)/fs; 
%end_time = 5;
end_time = find(y > 0.002, 1,'last')/fs; 
figure(15)
plot(start_time:1/fs:end_time,y(start_time*fs:end_time*fs));
title('Segmentirana govorna sekvenca'); 
xlabel('t[s]'); 
ylabel('y(t)');
xlim([start_time, end_time]);
 
%% Procena Tp preko autokorelacije

Xmax = max(abs(y));
Cl = 0.3*Xmax;

% Clipping
y_cl = zeros(length(y),1);
for i = 1:length(y)
    if y(i)>Cl
        y_cl(i) = 1;
    elseif y(i)<Cl
        y_cl(i) = -1;
    end
end
%y_cl(y>Cl) = 1;
%y_cl(y<-Cl) = -1;

n = 1000*8:6000*8;
figure;
    plot(t(n),y(n),'b',t(n),y_cl(n),'r');
    legend('Original signal', 'Clipped signal');
    title('Original vs clipped signal')
    xlabel('t[s]');

wl = fs*30e-3 + 1;
win = hamming(wl);
Ry = zeros(size(y_cl));

for i = (wl+1)/2: length(y_cl)-wl/2
    rng = (i-(wl-1)/2) : (i+wl/2); % centered
    rng0 = ((wl+1)/2-(wl-1)/2) : ((wl+1)/2+wl/2);
    Ry(i) = sum(y_cl(rng).*y_cl(rng0).*win.^2);
end
Ry = Ry/wl;

figure
    plot(t(n),y_cl(n),'b',t(n),Ry(n),'r');
    legend('Clipped signal', 'Autocorrelation function');
    title('Autocorrelation on clipped signal')
    xlabel('t[s]');

%% Konacna procena Tp

[~,indexes] = findpeaks(Ry);
indexes = indexes/fs;
Tp = [];
for i = 2:length(indexes)
    j = indexes(i)-indexes(i-1);
    if (j<1/100)&&(j>1/300)                
        Tp = [Tp j];
    end
end

Fp = 1/median(Tp);
disp(['Dobijena pitch frekvencija je: ' + string(Fp) + ' Hz']);
