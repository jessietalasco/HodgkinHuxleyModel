%% Simulation of linear RC circuit using forward eulers
% y 1-3 represent versions of the RC circuit sovled using different step
% sizes
% Initialize vectors to hold output voltage and add initial values
y1=zeros(1,100);
y1(1)=2;
y2=zeros(1,100);
y2(1)=2;
y3=zeros(1,100);
y3(1)=2;
% set up values for r and c
r=100*10^3;
c=50*10^(-6);
% set up input
x=zeros(1,100);
x(1,10:40)=2;
%set up different step sizes
h1=1/10*r*c;
h2=1/2*r*c;
h3=r*c;
% loop through input and calculate output using forward eulers method
tic
for i=2:100
    y1(i)=h1/c*x(i-1)+y1(i-1)*(1-h1/(r*c));
    y2(i)=h2/c*x(i-1)+y2(i-1)*(1-h2/(r*c));
    y3(i)=h3/c*x(i-1)+y3(i-1)*(1-h3/(r*c));
    
end
toc
% plot the output for each different step size
plot(y1)
hold on 
plot(y2)
hold on
plot(y3)
legend('h=tau/10','h=tau/2','tau')
xlabel('t')
ylabel('Vm')
title('Simulation of RC Circuit with Forward Eulers')
%% Simulation of linear RC circuit using backward eulers
% y1-3be represent the output of the RC circuit solved using different time
% steps
% initialize vectors to hold output and set initial values
y1be=zeros(1,100);
y1be(1)=2;
y2be=zeros(1,100);
y2be(1)=2;
y3be=zeros(1,100);
y3be(1)=2;
% initialize values of r and c
r=100*10^3;
c=50*10^(-6);
% set up input vector
x=zeros(1,100);
x(1,10:40)=2;
% iinitialize different time steps
h1=1/10*r*c;
h2=1/2*r*c;
h3=r*c;
% loop through the input and calculate the output using backward eulers for
% each step size.
tic
for i=2:100
    y1be(i)=r*h1/(r*c+h1)*x(i)+r*c/(r*c+h1)*y1be(i-1);
    y2be(i)=r*h2/(r*c+h2)*x(i)+r*c/(r*c+h2)*y2be(i-1);
    y3be(i)=r*h3/(r*c+h3)*x(i)+r*c/(r*c+h3)*y3be(i-1);
end
toc
% plot the output for each step size
plot(y1be)
hold on 
plot(y2be)
hold on
plot(y3be)
legend('h=tau/10','h=tau/2','tau')
xlabel('t')
ylabel('Vm')
title('Simulation of RC Circuit with Backward Eulers')
%% Plot comparison of forward and backward
% create new figure with subplots, plot each output with different time
% steps using each method with forward eulers on top and backward eulers on
% bottom, add legend, add axis labelsl, add titles.
figure
subplot(3,1,1)
plot(x)
xlabel('t')
ylabel('I')
title('Input')
subplot(3,1,2)
plot(y1)
hold on 
plot(y2)
hold on
plot(y3)
legend('h=tau/10','h=tau/2','tau')
xlabel('t')
ylabel('Vm')
title('Simulation of RC Circuit with Forward Eulers')
subplot(3,1,3)
plot(y1be)
hold on 
plot(y2be)
hold on
plot(y3be)
legend('h=tau/10','h=tau/2','tau')
xlabel('t')
ylabel('Vm')
title('Simulation of RC Circuit with Backward Eulers')
%% Using ODE45
% set up duration of simulation
sim_duration = 100; % sec
y_init = 0; % initial condition


options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')
tic % for measuring time it takes 
[time, y_sim] = ode45(@ydiff,[0 sim_duration],y_init,options); % pass in function describing derivative, time range, initial value, and options
toc
% create input vector for plotting
x=zeros(size(time));
for i=1:length(time)
    if time(i)>10 && time(i)<40
        x(i)=2;
    else 
        x(i)=0;
    end
end
% plot the output of ode45 and the input
figure
hold on
subplot(2,1,2)
plot(time, y_sim,'k','linewidth',2) % add RK4 result to plot
plot(time, y_sim,'k*','linewidth',2) % add RK4 result to plot
xlabel('t')
ylabel('voltage')
title('Using ODE45 to solve RC Circuit')
disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display
subplot(2,1,1)
plot(time,x)
xlabel('t')
ylabel('I')
title('Input')


function dydt = ydiff(time,y)
% diff EQ for ode45, for RC ckt problem
% Pass in the array y and the current time point, time

%initialize r and c
r=100*10^3;
c=50*10^(-6);
% set up the vector to hold the derivative
dydt = zeros(1,1);  % dimension the output variable; dydt = zeros(3,1)  for 3 diff EQs (column vector)
% create the input
if time >= 10 && time <= 40  % determine value of x, depending on time (secs)
    x = 2;
else
    x = 0;
end
dydt(1) = 1/c*x-1/(r*c)*y; % Differential equation describing the system
end
%% Cell Membrane with Constant Channel Conductances
%initialize simulation duration
sim_duration = 10; % sec
%set initial condition
y_init = 20.548; % initial condition


options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')

[time, y_sim] = ode45(@ydif,[0 sim_duration],y_init,options); % pass in function describing derivative, time range, initial value, and options
%set up input vector for plotting
x=zeros(size(time));
for i=1:length(time)
    if time(i)>2 && time(i)<8
        x(i)=2;
    else 
        x(i)=0;
    end
end

%plot output of ODE45 and input
figure
hold on
subplot(2,1,1)
plot(time,x)
xlabel('t')
ylabel('I')
title('Input')
subplot(2,1,2)
plot(time, y_sim,'k','linewidth',2) % add RK4 result to plot
plot(time, y_sim,'k*','linewidth',2) % add RK4 result to plot
title('Membrane Voltage')
disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display

function dydt=ydif(time,y)
% initialize vector to hold derivative
dydt=zeros(1,1);
% set nernst potentials, conductance, capacitance
Ena=50;
Ek=-77;
El=-54.4;
cm=1;
gna=120;
gk=36;
gl=0.3;
% create input
if time > 2 && time < 8
    x=2;
else 
    x=0;
end
%specify derivative
dydt(1)=(x-gna*(y-Ena)-gk*(y-Ek)-gl*(y-El))/cm;
end
%% Cell Membrane with Changing Channel Conductances
% set simulation duration
sim_duration = 100; % msec
%set all initial conditions
y_init = -64.9997; % initial condition
m0=0.0529625;
h0=0.59826;
n0=0.317764;
init_conds=[y_init;m0;h0;n0];

options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')

[time, y_sim] = ode15s(@ydiff2,[0 sim_duration],init_conds,options); % pass in function describing derivative, time range, initial value, and options
%initialize input for plotting 
x=zeros(size(time));
for i=1:length(time)
    if time(i) > 30 && time(i) < 33
        x(i)=20;
    elseif time(i)>38 && time(i)<41
        x(i)=20;
    elseif time(i)>48 && time(i)<51
        x(i)=20;
    elseif time(i)>55 && time(i)<58
        x(i)=20;
    elseif time(i)>65 && time(i)<68
        x(i)=20;
    elseif time(i)>71 && time(i)<74
        x(i)=20;
    else 
        x(i)=0;
    end

end

% plot input, output voltage, and activation/inactivation variables
figure
hold on
subplot(3,1,1)
plot(time,x) %plot x
title('Input Current')
xlabel('Time (ms)')
ylabel('Current(uA)')
subplot(3,1,2)
plot(time, y_sim(:,1),'k','linewidth',2) % add RK4 result to plot
plot(time, y_sim(:,1),'k*','linewidth',2) % add RK4 result to plot
title('Membrane Potential')
xlabel('Time (ms)')
ylabel('Voltage(mV)')
subplot(3,1,3)
plot(time,y_sim(:,2)) %plot m
hold on
plot(time,y_sim(:,3)) %plot h
hold on
plot(time,y_sim(:,4))% plot n
hold on
title('Channel Activation')
xlabel('Time (ms)')
ylabel('Probability')
legend('m','h','n')
disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display




function dydt=ydiff2(time,y)
%initialize vector to hold derivatives
dydt=zeros(4,1);
%set nernst potentials, max conductances, capacitance
Ena=50;
Ek=-77;
El=-54.4;
cm=1;
gna=120;
gk=36;
gl=0.3;
%set input 
if time > 30 && time < 33
    x=20;
elseif time>38 && time<41
    x=20;
elseif time>48 && time<51
    x=20;
elseif time>55 && time<58
    x=20;
elseif time>65 && time<68
    x=20;
elseif time>71 && time<74
    x=20;
else 
    x=0;
end

% set variable names to ys for ease of understanding
vm=y(1);
m=y(2);
h=y(3);
n=y(4);
% set alpha and beta equations for m
am=0.1*(vm+40)/(1-exp(-(vm+40)/10));
bm=4*exp(-(vm+65)/18);
% set alpha and beta equations for h
ah=0.07*exp(-(vm+65)/20);
bh=1/(1+exp(-(vm+35)/10));
% set alpha and beta equations for n
an=0.01*(vm+55)/(1-exp(-(vm+55)/10));
bn=0.125*exp(-(vm+65)/80);
% set derivative equations for voltage, m, h, and n
dydt(1)=(x-gna*m^3*h*(vm-Ena)-gk*n^4*(vm-Ek)-gl*(vm-El))/cm;
dydt(2)=am*(1-m)-bm*m;
dydt(3)=ah*(1-h)-bh*h;
dydt(4)=an*(1-n)-bn*n;
end
%% Simulating TTX
% Changed gna to 0 and set input back to step up to 20 from 30-70ms
% set dimulation duration
sim_duration = 100; % msec
%set initial conditions
y_init = -64.9997; % initial condition
m0=0.0529625;
h0=0.59826;
n0=0.317764;
init_conds=[y_init;m0;h0;n0];

options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')

[time, y_sim] = ode15s(@ydiff3,[0 sim_duration],init_conds,options); % pass in function describing derivative, time range, initial value, and options
%set input for plotting
x=zeros(size(time));
for i=1:length(time)
    if time(i)>30 && time(i)<70
        x(i)=20;
    else 
        x(i)=0;
    end
end

figure
hold on
subplot(3,1,1)
plot(time,x) % plot input current
title('Input Current')
xlabel('Time (ms)')
ylabel('Current(uA)')
subplot(3,1,2)
plot(time, y_sim(:,1),'k','linewidth',2) % add RK4 result to plot
plot(time, y_sim(:,1),'k*','linewidth',2) % add RK4 result to plot
title('Membrane Potential')
xlabel('Time (ms)')
ylabel('Voltage(mV)')
subplot(3,1,3)
plot(time,y_sim(:,2)) % plot m
hold on
plot(time,y_sim(:,3))% plot h
hold on 
plot(time,y_sim(:,4)) % plot n
hold on
title('Channel Activation')
xlabel('Time (ms)')
ylabel('Probability')
legend('m','h','n')
disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display




function dydt=ydiff3(time,y)
%initialize vector for derivatives
dydt=zeros(4,1);
% set nernst potentials, max conductances, and capacitance
Ena=50;
Ek=-77;
El=-54.4;
cm=1;
gna=120;
gk=36;
gl=0.3;
% set up input current
if time > 30 && time < 70
    x=20;
else 
    x=0;
end
% set variables to y values for ease of understanding
vm=y(1);
m=y(2);
h=y(3);
n=y(4);
% create alpha and beta equations for m
am=0.1*(vm+40)/(1-exp(-(vm+40)/10));
bm=4*exp(-(vm+65)/18);
% create alpha and beta equations for h
ah=0.07*exp(-(vm+65)/20);
bh=1/(1+exp(-(vm+35)/10));
% create alpha and beta equations for n
an=0.01*(vm+55)/(1-exp(-(vm+55)/10));
bn=0.125*exp(-(vm+65)/80);
% setting up derivative equations for solving
dydt(1)=(x-0*m^3*h*(vm-Ena)-gk*n^4*(vm-Ek)-gl*(vm-El))/cm;
dydt(2)=am*(1-m)-bm*m;
dydt(3)=ah*(1-h)-bh*h;
dydt(4)=an*(1-n)-bn*n;
end
%% Svirskis et al. model
% set simulation duration
sim_duration = 200; % msec
% specify initial conditions
y_init = -50; % initial condition
m0=0.0655;
h0=0.2352;
n0=0.0864;
k0=0.4005;
init_conds=[y_init;m0;h0;n0;k0];

options = odeset('RelTol',1e-6);  % Select Tolerance, to guide adaptive time step (try 1e-6; try 1e-2 - see 'chatter')

[time, y_sim] = ode15s(@ydiff4,[0 sim_duration],init_conds,options); % pass in function describing derivative, time range, initial value, and options
% create input for plotting
x=zeros(size(time));
for i=1:length(time)
    if time(i)>30 && time(i)<33
        x(i)=20;
    elseif time(i)>34 && time(i)<37
        x(i)=20;
    else 
        x(i)=0;
    end
end

figure
hold on
subplot(3,1,1)
plot(time,x) % plot input
title('Input Current')
xlabel('Time (ms)')
ylabel('Current(uA)')
subplot(3,1,2)
plot(time, y_sim(:,1),'k','linewidth',2) % add RK4 result to plot
plot(time, y_sim(:,1),'k*','linewidth',2) % add RK4 result to plot
title('Membrane Potential')
xlabel('Time (ms)')
ylabel('Voltage(mV)')
subplot(3,1,3)
plot(time,y_sim(:,2)) % plot m
hold on
plot(time,y_sim(:,3)) % plot h
hold on
plot(time,y_sim(:,4)) % plot n
hold on
plot(time,y_sim(:,5)) % plot activation of new channel
title('Channel Activation')
xlabel('Time (ms)')
ylabel('Probability')
legend('m','h','n','k')
disp(['# steps taken by ode45: ' num2str(length(time))]) % how many steps were used?
disp('  '); %Blank line to clarify display




function dydt=ydiff4(time,y)
% create vector to hold derivative values
dydt=zeros(5,1);
%specify nernst potentials, max conductances, and capacitance
Ena=50;
Ek=-90;
El=-54.4;
cm=1;
gna=0.2;
gk=0.01;
gkk=0.02;
gl=0.00333;
%set input
if time > 30 && time < 33
    x=20;
elseif time>34 && time<37
    x=20;
else 
    x=0;
end
% give ys clearer variable names
vm=y(1);
m=y(2);
h=y(3);
n=y(4);
k=y(5);

%set alpha and beta equations for m
am=4.2*exp(-0.0393*3.3*0.7*(-29.5-vm));
bm=4.2*exp(0.0393*3.3*(1-0.7)*(-29.5-vm));
%calculate m steady state
minf=am/(am+bm);
% calculate tau
taum=1/(am+bm);
taummin=0.05;
taum=max(taum,taummin); %ensure tau doesn't go below min

%set alpha and beta equations for h
ah=0.09*exp(-0.0393*-3*0.27*(-60-vm));
bh=0.09*exp(0.0393*-3*(1-0.27)*(-60-vm));
% calculate steady state value of h
hinf=ah/(ah+bh);
tauh=1/(ah+bh);
% calculate tau
tauhmin=0.25;
tauh=max(tauh,tauhmin); %ensure tau doesnt drop below min

%set alpha and beta equations for n
an=0.3*exp(-0.0393*3*0.8*(-30-vm));
bn=0.3*exp(0.0393*3*(1-0.8)*(-30-vm));
% calculate steady state value of n
ninf=an/(an+bn);
%set tau
taun=1/(an+bn);
taunmin=1;
taun=max(taun,taunmin); % ensure tau doesn't drop below min

%set alpha and beta equations for new channel
ak=0.2*exp(-0.0393*2.88*0.39*(-45-vm));
bk=0.17*exp(0.0393*2.88*(1-0.39)*(-45-vm));
%calculate steady state value
kinf=ak/(ak+bk);
%calculate tau (no minimum given)
tauk=1/(ak+bk);

%set derivative equations for voltage, m,h,n,and new channel
dydt(1)=(x-gna*m^3*h*(vm-Ena)-gk*n^4*(vm-Ek)-gkk*k*(vm-Ek)-gl*(vm-El))/cm;
dydt(2)=(minf-m)/taum;
dydt(3)=(hinf-h)/tauh;
dydt(4)=(ninf-n)/taun;
dydt(5)=(kinf-k)/tauk;
end