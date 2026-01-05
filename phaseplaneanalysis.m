%% Jessie Talasco Phase Plane Analyses
sim_duration = 100; % msec
y_init = -64.9997; % initial condition for voltage
w0=0.4017; % initial condition for w

init_conds=[y_init;w0]; %creating column vector of initial conditions

options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')

[time, y_sim] = ode15s(@ydiff2,[0 sim_duration],init_conds,options); % pass in function describing derivative, time range, initial value, and options
x=zeros(size(time));% setting up input current for graphing
for i=1:length(time)
    if time(i) > 25 && time(i) < 28
        x(i)=15;
    elseif time(i) > 29 && time(i) < 31
        x(i)=15;
    elseif time(i) > 32 && time(i) < 35
        x(i)=15;
    elseif time(i) > 36 && time(i) < 39
        x(i)=15;
    elseif time(i) > 40 && time(i) < 43
        x(i)=15;
    else
        x(i)=0;
    end
end

figure
hold on

subplot(3,1,1)
plot(time,x) %graph input current
title('Input Current')
xlabel('Time (ms)')
ylabel('Current(uA)')

subplot(3,1,2) 
plot(time, y_sim(:,1),'k','linewidth',2) % Graph voltage vs time
plot(time, y_sim(:,1),'k*','linewidth',2) 
title('Membrane Potential')
xlabel('Time (ms)')
ylabel('Voltage(mV)')

subplot(3,1,3)
plot( y_sim(:,1), y_sim(:,2)) % graph w vs V (Phase Plot)
title('Phase Plane')
xlabel('Voltage (mV)')
ylabel('w')

disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display




function dydt=ydiff2(time,y)
dydt=zeros(2,1); %initiate vector for holding derivative values
% set values of nernst potentials, membrane capacitance, and max
% conductances
Ena=50;
Ek=-77;
El=-54.4;
cm=1;
gna=120;
gk=36;
gl=0.3;
% set up stimulus
if time > 25 && time < 28
    x=15;
elseif time > 29 && time < 31
    x=15;
elseif time > 32 && time < 35
    x=15;
elseif time > 36 && time < 39
    x=15;
elseif time > 40 && time < 43
    x=15;
else
    x=0;
end
vm=y(1); %name v for ease of reading
w=y(2); %name w for ease of reading
s=(1-0.59826)/0.317764; %calculate s based off of previous initial conditions for n and h
% calculate alphas and betas from given formulas
am=0.1*(vm+40)/(1-exp(-(vm+40)/10)); 
bm=4*exp(-(vm+65)/18);
ah=0.07*exp(-(vm+65)/20);
bh=1/(1+exp(-(vm+35)/10));
an=0.01*(vm+55)/(1-exp(-(vm+55)/10));
bn=0.125*exp(-(vm+65)/80);
%calculate the steady state values using alphas and betas
ninf=an/(an+bn);
hinf=ah/(ah+bh);
minf=am/(am+bm);
winf=s*(ninf+s*(1-hinf))/(1+s^2);
%calculate the tau for w using the formula giben in Rinzel
tauw=5*exp(-(vm+10)^2/55^2)+1;


dydt(1)=(x-gna*minf^3*(1-w)*(vm-Ena)-gk*(w/s)^4*(vm-Ek)-gl*(vm-El))/cm; %diff eq for voltage
dydt(2)=3.8*(winf-w)/tauw; %diff eq for w
end