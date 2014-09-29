% Hodgkin Huxley Model- Phase I

%%% Constants

gK=36; %Potassium channel conductance (mS/cm^2)
gNA=120; %Sodium channel conductance (mS/cm^2)
gL=0.3; %Leakage channel conductance (mS/cm^2)
Ek=-12; %volt across K channel (mV)
Ena=115; %Volt across Na channel (mV)
El=10.6; %Volt across leakage channel (mV)
V_rest=0; %resting membrane voltage (mV) (we will later subtract by 70 to 
%receive accurate values. Code does not run properly with V_rest= -70
Cm=1; %Membrane capacitance (microFarads/cm^2)


%%% Initial time/other initial conditions

%Total simulation time: 100ms
t=0; % initial time (ms)
ss=1; %Step size
Vm= V_rest; %initial voltage (mV)

% Initialize Gating variables:
Am=0.1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm=4*exp(-Vm/18);
An=.01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn=.125*exp(-Vm/80);
Ah=.07*exp(-Vm/20);
Bh=1/(exp((30-Vm)/10)+1);
m=Am/(Am+Bm);
n=An/(An+Bn);
h=Ah/(Ah+Bh);

% Currents:
I=0; %Steady state neuron will have no current
Il=gL*(Vm-El); %Leakage current
Ik=n^4*gK*(Vm-Ek); %Current through K channels
Ina=m^3*gNA*h*(Vm-Ena); %Current through Na channels
I_ion=I-Ina-Ik-Il;

% Derivatives:
dVm_dt=I_ion/Cm; %change in volt/change in time
Vm= Vm+ss*dVm_dt; %updated voltage

dm_dt=Am*(1-m)-Bm*m; 
dn_dt=An*(1-n)-Bn*n;
dh_dt=Ah*(1-h)-Bh*h;

%Update m, n, h
m= m+ss*dm_dt;
n= n+ss*dn_dt;
h= h+ss*dh_dt;

%Vectors (used later to plot graph)
VecX= (0:100); %Time vector
VecY= (Vm); %Voltage vector
VecgK= (gK); %Potassium conductance vector
VecgNA= (gNA); %Sodium conductance vector

%Update time
t= t+ss;

%Loop (# of elements~ time)
while t <= 100
    
    % Gating variables:
    Am=0.1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm=4.*exp(-Vm/18);
    An=.01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn=.125.*exp(-Vm/80);
    Ah=.07*exp(-Vm/20);
    Bh=1/(exp((30-Vm)/10)+1);

    m=Am/(Am+Bm);
    n=An/(An+Bn);
    h=Ah/(Ah+Bh);

    % Currents:
    I=0;
    Il=gL.*(Vm-El);
    Ik=(n.^4).*gK.*(Vm-Ek);
    Ina=(m.^3).*gNA.*h.*(Vm-Ena);
    I_ion=I-Ina-Ik-Il;

    % Derivatives for updates:
    dVm_dt=I_ion/Cm;
    dm_dt=Am*(1-m)-Bm*m;
    dn_dt=An*(1-n)-Bn*n;
    dh_dt=Ah*(1-h)-Bh*h;
    
    % Update m,n,h
    m=m+(ss*dm_dt);
    n=n+(ss*dn_dt);
    h=h+(ss*dh_dt);
    
    % Euler's method-update voltage
    Vm= Vm+ss*dVm_dt;
    
    % Update time
    t= t+ss; 
    
    % Adding new values to voltage & conductance vectors
    VecY= [VecY Vm]; %Voltage vector
    VecgNA= [VecgNA gNA]; %Sodium conductance vector
    VecgK= [VecgK gK]; %Potassium conductance vector
    
end

VecY= VecY-70; %This is to correct for starting V_rest at 0

%Plotting
plot(VecX, VecY, 'c-')
title('Action Potential: Change of Voltage with Time')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis([0,100,-100,150])
