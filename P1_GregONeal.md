BMED-4784-Project
=================
%Initialize time
%%%Total simulation time: 100ms
t=1;

%Constants
gK=36;
gNA=120;
gL=0.3;
Ek=-12;
Ena=115;
El=10.6;
V_rest=-70;
Cm=1;

% Initial conditions
Vm=V_rest;
VecX= (0:100);
VecY= [-70];


%Loop (# of elements symbolize time)
while t <= 100
    
    % Gating variables:
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
    Il=gL*(Vm-El);
    Ik=n^4*gK*(Vm-Ek);
    Ina=m^3*gNA*h*(Vm-Ena);
    I_ion=Ina+Ik+Il;

    % Derivatives:
    dVm_dt=I_ion/Cm;
    dm_dt=Am*(1-m)-Bm*m;
    dn_dt=An*(1-n)-Bn*n;
    dh_dt=Ah*(1-h)-Bh*h;
    

    ss= 1; %step size
    
    % Update m,n,h
    m=m+(ss*dm_dt);
    n=n+(ss*dn_dt);
    h=h+(ss*dh_dt);
    
    %Euler's method-update voltage and time
    Vm= Vm+ss*dVm_dt;
    t= t+1; 
    
    %Vector of 101 voltage points (changing over time)
    VecY= [VecY Vm];
    
end


%Plotting
plot(VecX, VecY, 'c-')

title('Action Potential: Change of Voltage with Time')
xlabel('Time (ms)')
ylabel('Voltage (mV)')
