 clc
clear all
close all
%This code solves the microkinetic model below giving the kinetic current
%--------simulation of a catalytic reaction
% P+e-->Q
%Q-->B
%B+e-->C
%C+e-->D
%D+e-->P


syms y1 y2 y3 y4 y5 En1 delK GO2 KOH k1 k2 k3 k4 k5 KO2 Kads E gamma

%----Parameters--------
 Z=load(''); %Load data in mA/cm^2

F=96485; %Faraday's constant
R=8.31; %Gas constant 
Tem=295.15; %Temperature in K
b=R.*Tem./F; %Thermal voltage in V
Cb=1.2E-6; %1.2E-6 Oxygen concentration in bulk, mol/cm^3
DO2=2E-5; %2E-5 diffusion coefficient in 0.1 M KOH cm^2/s
%-----Experimental Paramteres----------------
L=0.561; %Diameter of electrode cma
A=0.2475;

alpha1=0.5;
alpha3=0.5;
alpha4=0.5;
alpha5=0.5;
k_2=k2./Kads;
En3=4.92-delK-(R.*Tem.*log(Kads)./F)-En1; 
En4=delK+2*KOH-En1;
En5=En1-2*KOH;
%----------Dimensional analysis-------------------------------
ts=L.^2./DO2;
%---Dimensional analysis for SS-----------------

E1=En1./b;
E3=En3./b;
E4=En4./b;
E5=En5./b;
V=E./b;
%---Define rate constants-----------------------
xsi1=(V-E1); xsi3=(V-E3); xsi4=(V-E4);xsi5=(V-E5); 
K1=k1.*exp(-xsi1*alpha1)*ts;
K_1=k1.*exp(xsi1*(1-alpha1))*ts;
K2=k2*ts;
K_2=k_2.*ts;
K3=k3.*exp(-xsi3*(alpha3)).*ts;
K_3=k3.*exp(xsi3*(1-alpha3)).*ts;
K4=k4.*ts.*exp(-xsi4*alpha4);
K_4=k4.*ts.*exp(xsi4*(1-alpha4));
K5=k5*ts.*exp(-xsi5*(alpha5));
K_5=k5.*ts.*exp(xsi5*(1-alpha5));
%----Define the set of equations---------------
r1=K1.*y1-K_1.*y2;
r2=K2.*y2-K_2.*y3;
r3=K3.*y3-K_3.*y4;
r4=K4.*y4-K_4.*y5;
r5=K5.*y5-K_5.*y1;

%----Solve set of equations using matrix solver--------------------
vars=[y1, y2, y3, y4, y5];

eq=[r5-r1==0, r1-r2==0, r2-r3==0, r3-r4==0,r4-r5==0,(y1+y2+y3+y4+y5)==1];
[M,S] = equationsToMatrix(eq,vars);

X=M\S;%Find surface coverage



%Insert surface coverage
R1=K1.*X(1)-K_1.*X(2);
R3=K3.*X(3)-K_3.*X(4);
R4=K4.*X(4)-K_4.*X(5);
R5=K5.*X(5)-K_5.*X(1);

Jk=-(R1+R3+R4+R5);%Find kinetic current
jk=matlabFunction(Jk); %convert to function

%----Insert values into the function
En1=0.818;
delK=3.03;
GO2=4.92;
KOH=-0.28;
k1=60;
k2=2000;
k3=3;
k4=50;
k5=50;
KO2=4.6058;
Kads=0.6;
E=linspace(1,0.55,1000);
gamma=1.13e-10;
J_Norm=gamma.*F./ts;

numjk=jk(E,En1,KOH,Kads,delK,k1,k2,k3,k4,k5).*J_Norm;
toc
plot(E,1000*numjk)
hold on
plot(Z(:,1),Z(:,2),'o');
