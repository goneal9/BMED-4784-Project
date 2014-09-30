% Hodgkin Huxley Model- Phase I

%%% Constants
maxCondK=36; %Potassium channel conductance (mS/cm^2)
maxCondNA=120; %Sodium channel conductance (mS/cm^2)
maxCondL=0.3; %Leakage channel conductance (mS/cm^2)
Ek=-12; %volt across K channel (mV)
Ena=115; %Volt across Na channel (mV)
El=10.6; %Volt across leakage channel (mV)
V_rest=0; %resting membrane voltage (mV) (we will later subtract by 70 to 
%receive accurate values. Code does not run properly with V_rest= -70
Cm=1; %Membrane capacitance (microFarads/cm^2)

%%% Initial time/other initial conditions
total_t=100; %total time= 100ms
t=0; % initial time (ms)
ss=.01; %Step size= .01
Vm= V_rest; %initial voltage (mV)

% Initialize Gating variables:
Am=0.1*((25-Vm)/(exp((25-Vm)/10)-1));
Bm=4*exp(-Vm/18);
An=.01*((10-Vm)/(exp((10-Vm)/10)-1));
Bn=.125*exp(-Vm/80);
Ah=.07*exp(-Vm/20);
Bh=1/(exp((30-Vm)/10)+1);

%Initial m,n,h (actvation probabilities)
m=Am/(Am+Bm); %probability that particle m is "open
n=An/(An+Bn); %probability that particle is "open"
h=Ah/(Ah+Bh); %probability that particle h is "open"

gNA= m.^3*h*maxCondNA; %Conductance for NA
gK= n.^4*maxCondK; %Conductance for K

%Vectors (used later to plot graph)
VecX= [0:.01:100]; %Time vector
VecY= (Vm); %Voltage vector
VecgK= (gK); %Potassium conductance vector
VecgNA= (gNA); %Sodium conductance vector


%Loop (# of elements~ time)
while t <= 100
    
    % Gating variables update:
    Am=0.1.*((25-Vm)/(exp((25-Vm)/10)-1));
    Bm=4.*exp(-Vm/18);
    An=.01*((10-Vm)/(exp((10-Vm)/10)-1));
    Bn=.125.*exp(-Vm/80);
    Ah=.07*exp(-Vm/20);
    Bh=1/(exp((30-Vm)/10)+1);
    
    % Currents:
    I=0; %steady state current=0
    Il=maxCondL.*(Vm-El); %leakage current (All other ions besides Na+ & K+)
    Ik=(n.^4).*maxCondK.*(Vm-Ek); %current across K channel
    Ina=(m.^3).*maxCondNA.*h.*(Vm-Ena); %Current across sodium channel
    I_ion=Ina+Ik+Il; %I_ion is the ionic current 
    Ic=I-I_ion;%Ic is the capacitative current
    
    % Update voltage derivative
    dVm_dt=Ic/Cm; 
    % Use Euler's method to update voltage
    Vm= Vm+ss*dVm_dt;
    
    % Change in activation probabilities
    dm_dt=Am*(1-m)-Bm*m; %movement of particle m from closed to open
    dn_dt=An*(1-n)-Bn*n; %movement of particle from closed to open
    dh_dt=Ah*(1-h)-Bh*h; %movement of particle h from closed to open
    
    %Update m, n, h
    m= m+ss*dm_dt;
    n= n+ss*dn_dt;
    h= h+ss*dh_dt;
    
    % Adding new values to voltage & conductance vectors
    VecY= [VecY Vm]; %Voltage vector
    
    %Updating conductance values:
    gNA= m.^3*h*maxCondNA; %Conductance for NA
    gK= n.^4*maxCondK; %Conductance for K
    VecgNA= [VecgNA gNA]; %Sodium conductance vector
    VecgK= [VecgK gK]; %Potassium conductance vector
    
    % Update time
    t= t+.01;  
end

VecY= VecY-70; %This is to correct for starting V_rest at 0

%Plotting 

%Voltage as a function of time
figure
Volt_time= plot(VecX, VecY, 'k-');
title('Membrane Potential: Change of Voltage over Time for Steady State Neuron')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
axis([0,100,-100,150])

%Sodium/Potassium conductance
figure
cNA= plot(VecX, VecgNA);
hold on
cK= plot(VecX, VecgK, 'r');
title('Channel Conductances of Na/K for Steady State Neuron')
legend([cNA,cK],'Conductance for Na+','Conductance for K+')
axis([0,100,0,5])
ylabel('conductance (mS/cm^2)')
xlabel('Time(ms)')
