%% Obrada i prepoznavanje govora zadatak 3


clear ;
close all; 
clc;

%% Skriveni Markovljev model - problem 1

a = 0.05;
b = 0.1;

A = [1-3*a a 2*a; b 1-2*b b; 0.1 0.1 0.8];
B = [5/8 2/8 1/8; 2/13 7/13 4/13; 1/10 3/10 6/10];
P = [1/3 1/3 1/3];

% Opservacije
N = 100;
[s, v] = generisiOpservacije(N,P,A,B);

% Problem1
P1 = forwardBackward(v,P,A,B);
disp(['Verovatnoca tacne rekonstrukcije stanja u sekvenci: ' + string(P1)])

% Problem2
[x, P2] = viterby(v,P,A,B);

tacno = sum(s==x);
disp(['Tacno pogodjenih stanja: ' + string(tacno) + '/' + string(N)])
%% Skriveni Markovljev model - problem 2

a = 0.2;
b = 0.33;

A = [1-3*a a 2*a; b 1-2*b b; 0.1 0.1 0.8];
B = [5/8 2/8 1/8; 2/13 7/13 4/13; 1/10 3/10 6/10];
P = [1/3 1/3 1/3];

% Opservacije
N = 100;
[s, v] = generisiOpservacije(N,P,A,B);

% Problem1
P1 = forwardBackward(v,P,A,B);
disp(['Verovatnoca tacne rekonstrukcije stanja u sekvenci: ' + string(P1)])

% Problem2
[x, P2] = viterby(v,P,A,B);

tacno = sum(s==x);
disp(['Tacno pogodjenih stanja: ' + string(tacno) + '/' + string(N)])

