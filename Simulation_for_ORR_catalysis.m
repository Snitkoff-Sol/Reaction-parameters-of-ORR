 clc
clear all
close all
%Simulation used to model the ORR electrocatalysis
%in Snitkoff-Sol et. al. 2023 Nature Catalysis


%----Constant Parameters--------
F=96485; %Faraday's constant
R=8.31; %Gas constant 
Tem=295.15; %Temperature in K
b=R.*Tem./F; %Thermal voltage in V
Cb=1.2E-6; %Oxygen concentration in bulk, mol/cm^3
DO2=2.0E-5; %2E-5 diffusion coefficient in 0.1 M KOH cm^2/s
%-----Experimental Paramteres----------------
Amp=150/1000; %Amplitude in V
v=-4.883/1000; %Scan rate V/s
f=4.883; %Frequency in Hz
Ei=0.95; %Initial potential in V
Ef=0.55; %End potential in V
t_eCA=linspace(0,10,2^14); %Time for the CA
L=0.5; %Radious of electrode cma
A=0.196;%Area of electrode
Ru=90; %Uncompensated resistance in Ohm
Cdl=0.16E-3; %Double layer capacitance F/cm^2
%-----Define Exponential Spatial Grid---------
N=64;%Number of grids
beta=0.01;%Control parameter - defines the "speed" of change
dx=zeros(N,1);
M=1E-8;%Control parameter - defines the smallest distance
for l=1:N
dx(l)=M*exp(beta*(l-1));
M=dx(l);
end

%---------------Call Experimental Data---------------
Exp=load("");%Load data E|i|t
[C,ia,ic]=unique(Exp(:,3)-Exp(1,3)); %Find unique points
Np=14; %Number of data points
t_e=linspace(0,C(end),2^Np);%Time stamps of measurement
ExpInt=interp1(C,Exp(ia,2),t_e'); %Interpolate to fill in missing points
df=length(t_e)./t_e(end); 
Y=zeros(length(t_e),4);
band=1;
%Extract harmonics
for i=1:6
H=i+1;
wind=[H*f-band/2,H*f+band/2];
YHarC=bandpass(ExpInt(:),wind,df,'ImpulseResponse','iir','Steepness',0.8);
Y(:,i)=abs(hilbert(YHarC));
end
%--------------Kinetic and Thermodynamic Parameters-------------
%All parameters are given in Table 2 in the main text
En1=0.818;
GO2=4.92;
delK=3.03;
KOH=-0.28;
k1=60;
k2=2000;
Kads=0.6;
k3=3;
k4=30;
k5=45;

k_2=k2./Kads;
En3=GO2-delK-(R.*Tem.*log(Kads)./F)-En1; 
En4=delK+2*KOH-En1;
En5=En1-2*KOH;
gamma=2.2e-10;

%-----Dimensionless Parameters--------------

ts=L.^2./DO2;
Time=t_e.*f; 
TimeCA=t_eCA.*f;
Xsii=Ei./b;
K1=k1./f;
K2=k2./f;
K_2=k_2./f;
K3=k3./f;
K4=k4./f;
K5=k5./f;
Am=Amp/b;
Ni=v./(f.*b);
E1=En1./b;
E3=En3./b;
E4=En4./b;
E5=En5./b;
J_Norm=f*gamma*F;
tao=ts*f;
RC=Ru.*Cdl.*A.*f;
eps=gamma*F./(b.*Cdl);
mu=DO2.*Cb./(L*f*gamma);
%----------Find Initial Conditions------------------
C0=zeros(N,1);
C0(:)=1;
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
tic
[t1,y1]=ode15s(@(t,y)ORR_Cat(t,y,K1,K2,K_2,K3,K4,K5,Ni,0,Xsii,E1,E3,E4,E5,dx,N,tao,mu,eps,RC),TimeCA,[1,0,0,0,0,Xsii, C0'],opts);
toc

%--------Run Full Simulation-----------------------
Initi=[y1(end,1),y1(end,2),y1(end,3),y1(end,4),y1(end,5),Xsii,y1(end,7:end)];
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
tic
[t,y]=ode15s(@(t,y)ORR_Cat(t,y,K1,K2,K_2,K3,K4,K5,Ni,Am,Xsii,E1,E3,E4,E5,dx,N,tao,mu,eps,RC),Time,[Initi],opts);
toc
E=Xsii+Ni.*t+Am.*sin(2*pi*t);
%------Simulated Total Current----------------------------
J(:)=((E(:)-y(:,6)).*b./Ru).*1000;
figure(1)
plot(E.*b,J),title('Total Current'), xlabel('E [V vs  RHE]'), ylabel('i [mA]')
%--------Extract Harmonics--------------------
tdc=Ei+t_e.*v;
I_sim=zeros(length(t_e),6);
for i=1:6
H=i+1;
wind=[H*f-band/2,H*f+band/2];
YHarC=bandpass(J,wind,df,'ImpulseResponse','iir','Steepness',0.6);
I_sim(:,i)=abs(hilbert(YHarC));
end
%-------------Plot Simulated and Experimental Results--------
for j=1:6
    figure(j)
    k=num2str(j+1);
    plot(tdc,I_sim(:,j),tdc,Y(:,j)),title([k '^{th} Harmonic']), legend('Simulation','Experiment'), xlabel('Cell Voltage/V'), ylabel('i/mA')
end


%-------- Numerical solution to the diffusion problem and microkinetic model

function dydt = ORR_Cat(t,y,K1,K2,K_2,K3,K4,K5,Ni,Am,Xsii,E1,E3,E4,E5,dx,N,tao,mu,eps,RC)


V=Xsii+Ni.*t+Am.*sin(2*pi*t);

C=zeros(N,1);
Ct=zeros(N,1);
reactions=6; %Number of non-diffusional equations to solve

for i=1:N
    C(i)=y(i+reactions);
end

C(N)=1; 
xsi1=(y(6)-E1)./(2); xsi3=(y(6)-E3)./(2); xsi4=(y(6)-E4)./(2); xsi5=(y(6)-E5)./(2);
y(1)=1-(y(2)+y(3)+y(4)+y(5));
r1=K1*(y(1).*exp(-xsi1)-y(2).*exp(xsi1));
r2=K2.*y(2).*C(1)-K_2.*y(3);
r3=K3*(y(3).*exp(-xsi3)-y(4).*exp(xsi3));
r4=K4*(y(4).*exp(-xsi4)-y(5).*exp(xsi4));
r5=K5*(y(5).*exp(-xsi5)-y(1).*exp(xsi5));
rE=(V-y(6))./RC+eps.*(r1+r3+r4+r5);

for i=1:N-1
    if i==1
        C_0=C(2)-2*dx(1)*r2./mu;
        Ct(1)=((C(2)-C(1))./dx(2)+(C_0-C(1))./dx(1))./(1./2*(dx(1)+dx(2)))./tao;
    elseif i==N-1
        Ct(i)=((1-C(N-1))./(dx(N-1))+(C(N-2)-C(N-1))./dx(N-2))./(1./2*(dx(N-1)+dx(N-2)))./tao;
    else
        Ct(i)=((C(i+1)-C(i))./dx(i+1)+(C(i-1)-C(i))./dx(i))./(1./2.*(dx(i+1)+dx(i)))./tao;
    end
end


dydt=[r5-r1;r1-r2;r2-r3;r3-r4;r4-r5;rE;Ct];
end

