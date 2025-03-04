%%
% 1D Levee Building Model
% Advection settling of suspended sediment
% Jeongyeon Han and Wonsuck Kim
% July, 2020
%% Test 2
% Water Level Case 2 (same water depth HI at the floodplain)
% Entrainment O
%% 
clear all;
%%
g = 9.81; %acceleration of gravity
nu = 1*10^(-6); %kinematic viscosity
ro = 1000; %water density
ros = 2650; %sediment density
R = (ros/ro) -1; %submerged specific gravity
%% Model Parameters
M = 20; %number of nodes
L = 200; %levee length [m]
%% Overflow Properties
Hf = 4; %initial flow depth [m]
Uf = 0.1; %initial flow velocity [m/s]
qw = Hf*Uf; % water discharge; always constant [m2/s]
%% Assign Matrix of Parameters
dx = L/M; %step length
dt = 1; %time variation [s]
Niterations = 700000; %number of total iterations 

qs = zeros(M,7); %sediment discharge
qst = zeros(M,1); %total sediment discharge of each node
cibar = zeros(M,7); % sediment concentration
qsoi = zeros(7,1); %initial sediemnt discharge of each sediment grainsize range

Pi = zeros(M,7); %fractional distribution of ith grain size at each node in the flood flow
Pdi = zeros(M,7); %fractional distribution of ith grain size at each node in the deposit
f = zeros(M,7); %cumulative distribution of ith suspended grain size at each node
fd = zeros(M,7); %cumulative distribution of ith deposit grain size at each node
Psibar = zeros(M,1); %mean grainsize in psi unit at each node
D50 = zeros(M,1); % 50% grain size at each node
D90 = zeros(M,1); % 90% grain size at each node

Ds = zeros(M,7); %deposit rate for each grain size on one iteration
Dst = zeros (M,1); %total rate at each node on one iteration
T_Ds = zeros(M,7); %total deposition rate for each grain size at each node
T_Dst= zeros(M,1); %total deposition rate at each node

% Time variation parameter
jprint = zeros(1,1); 
jj = 0;

%% Assign Initial Slope
eta = zeros (M,1); %bed elevation
SI = 0.0001; %initial slope
etaI = L*SI; %inital elevation at levee front

%% Setup Initial Bed Surface (substrate)
for i=1:M
    eta(i) = etaI - SI*dx*i; %initial elevation at each node
end

%% Entrainment Relationship by Garcia and Parker (1991)
Aa = 1.3e-7; %Entrainment relationship, garcia and Parker (1991), A
Slope = zeros(M,1); %slope of levee
Zu = zeros(M,7);    %similarity collapse of each grain size at each node
Eui = zeros(M,7);   % entrainment from the bed of each grain size 
Hs = zeros(M,1);    % water depth associated with skin friction 
Cfn1o2 = zeros(M,1);    % Cf^-1/2
Ustars = zeros(M,1);    % U*s shear velocity of skin friction
Ei = zeros(M,7);    % total entrainment of each grain size at each node
T_Es = zeros(M,7);
Et = zeros(M,1);
%%
% Grain size 
psi = [-6 -5 -4 -3 -2 -1 0]; %sediment grain size
D = 2.^(psi); %sediment grain size in mm
D = D/1000; %sediment grain size in m
psibar = mean(psi);     %-3; %mean psi
sigma = 0.8;    %standard deviation for normal distribution
p =(1/(sqrt(2*pi)*sigma))*exp(-(((psi-psibar)./sigma).^2)/2); %gaussian grain size distribution
F = cumsum(p);
% Settling velocity of each grain size using Ferguson and Church (2004)
% The constants, C1 = 18 and C2 = 1 for natural grains
ws=(R*g.*(D.^2))./(18*nu+(0.75*R*g.*D.^3).^0.5);
Repi = ((R*g*D).^(0.5)).*D/nu;  % particle Reynolds Number

%% Rouse-Vanoni Equilibrium Suspended sediment profile
ca =(0.00016).*p; %near bed concentration
Hc = 8; % flood channel depth
Uc = 1.5; % channel velocity
Sc = 0.0001; % channel bed slope
ustar_c = sqrt(g*Hc*Sc); %channel shear velocity
P = ws./(0.4*ustar_c); %Rouse number: the exponent of Rouse equation
cibar0H = zeros(7,1);
cibarH = zeros(7,1);
T_qso = zeros(7,1);

%% Calculate Overflow Concentration
for k=1:7
fun = @(z) ca(k)*(((Hc./z - 1)./19).^P(k)); % Rouse equation
cibar0H(k) = integral(fun,(Hc-Hf), Hc); %input concentration
cibarH(k) = integral(fun,0.05*Hc,Hc);
T_qso(k) = Uc*(cibarH(k)+ca(k));
qsoi(k) = Uc*cibar0H(k); %initial sediment discharge
end
T_Q=sum(T_qso); %total sediment discharge in channel
T_qs=sum(qsoi);

%% Initial Concentration of Suspended Sediment
alpha = 1; % constant for initial suspended sediment concentration, 0 for fresh water
for i=1:M
    cibar(i,:) = (alpha.*qsoi)./(Uf*Hf);
end

%% Update Suspended Sediment Concentration and Levee Elevation
for j=1:1:Niterations
	for i=1:1:M
        for k=1:1:7 %7 grain size
                if i==1
                    cibar(i,k) = cibar(i,k)+ (((qsoi(k)-qs(i,k))/dx)+ws(k)*(Ei(i,k)-cibar(i,k)))*dt/Hf; %sediment concentration at first node using the ghost node
                else
                    cibar(i,k) = cibar(i,k)+ (((qs(i-1,k)-qs(i,k))/dx)+ws(k)*(Ei(i,k)-cibar(i,k)))*dt/Hf;  %sediment concentration at other nodes
                end
                qs(i,k) = cibar(i,k)*Uf*Hf; %sediment discharge at each node
                Ds(i,k) = ws(k)*(cibar(i,k)-Ei(i,k)); %deposit rate  for each grain size at each node on one iteration 

        end %for k
                qst(i) = sum(qs(i,:)); %total sediment discharge at each node 
                Dst(i) = sum(Ds(i,:)); %total deposit rate at each node on one iteration
                eta(i) = eta(i) + Dst(i)*dt; %update sediment surface elevation (Exner equation)
                
                %calculate slope
                if i == 1
                    Slope(i) = (eta(i)-eta(i+1))/dx;
                else
                    Slope(i) = (eta(i-1)-eta(i))/dx;
                end

                Pi(i,:) = qs(i,:)./qst(i);  % fraction function of each suspended grain size at each node
                %change deposition rate with entrainment
                Dstt = 0;
                for k=1:1:7
                   if Ds(i,k)>0
                       Dstt = Dstt + Ds(i,k);
                   end
                end   
                
                % using the total depositional rate compared to individual
                % depositional rate to calculate fraction of ith grainsize
                % range in the deposit
                for k=1:1:7
                    if Ds(i,k)>0
                        Pdi(i,k) = Ds(i,k)/Dstt;
                    else
                        Pdi(i,k) = 0;
                        Ds(i,k) = 0;
                    end
                    T_Ds(i,k) = T_Ds(i,k) + Ds(i,k); %total deposital rate at each node %일정 간격마다 추적하여 쌓이는 과정 plot하여 살피기
                 end

                T_Dst(i) = sum(T_Ds(i,:)); %total deposital rate at each node 
                f = cumsum(Pi')';   % cumulative distribution at each node for suspended sediment
                fd = cumsum(Pdi')'; % cumulative distribution at each node for deposit
                f(f>1)=1;   % some errors, sum larger than 1 changes to 1
                fd(fd>1)=1;   % some errors, sum larger than 1 changes to 1

                %calculating D50 
                xtx = 0.5; % grain size such that xtx % of the material
                u = 0;
                for y = 2:1:7
                    if fd(i,y) >= xtx && fd(i,y-1) < xtx
                      u = y; 
                    end
                end
                psix = psi(u-1) + ((psi(u)-psi(u-1))/(fd(i,u)-fd(i,u-1)))*(xtx - fd(i,u-1));  % psi values of x% of material
                D50(i,1) = 2^(psix); %grain size values of x% of material
                D50(i,1) = D50(i,1)/1000; %m

                %calculating D90
                xtx = 0.9; % grain size such that xtx % of the material
                u = 0;
                for y = 2:1:7
                    if fd(i,y) >= xtx && fd(i,y-1) < xtx
                      u = y; 
                    end
                end
                psix = psi(u-1) + ((psi(u)-psi(u-1))/(fd(i,u)-fd(i,u-1)))*(xtx - fd(i,u-1));  % psi values of x% of material
                D90(i,1) = 2^(psix); %grain size values of x% of material
                D90(i,1) = D90(i,1)/1000; %m

                %input Entrainment relationship by garcia and Parker (1991)
                Psibar(i) = sum(Pdi(i,:)*psi(:)); %mean psi values
                sigma2 = sum((psi-Psibar(i)).^2.*Pdi(i,:));
                sigma = (sigma2)^0.5; %standard deviation
                lamda_m = 1-0.298*sigma;    % straining parameter
                ks = 2*D90(i,1); %roughness height over a flat bed 
                Hs(i,1) = ((Uf*ks^(1/6))/(8.32*sqrt(g*Slope(i,1))))^(3/2);%mean depth associated with skin friction
                Cfn1o2(i,1) = 1/0.4*log(11*Hs(i,1)/ks); %friction coefficient
                Ustars(i,1) = Uf/Cfn1o2(i,1); %shear velocity
                Zu(i,:) = lamda_m * Ustars(i,1)./ws'.*((Repi').^0.6).*((D')./D50(i)).^0.2; 
                Eui(i,:) = (Aa.*Zu(i,:).^5)./(1+Aa/0.3.*(Zu(i,:).^5)); %dimensionless entrainment
                Ei(i,:) = Eui(i,:).*Pdi(i,:); %dimensionless entrainment rate per unit area
                
                for k=1:1:7
                    T_Es(i,k) = T_Es(i,k) + ws(k)*Ei(i,k); %total deposital rate at each node 
                 end

                Eit(i) = sum(T_Es(i,:)); %total deposital rate at each node 
                
	end % for M
        
        T_Pi = T_Ds./T_Dst;
        T_f = cumsum(T_Pi')'; 
        
         if eta(1) >= 2
            fprintf("\n avulsion time : %g sec \n", j*dt);
            break;
         end
        
      %print every 500 iteration
      if mod(j,500) 
        continue
      end
   jj = jj+1;   
   jprint(jj,1) = j;
    %upstream time variation (M=3)
   up_D50(jj,1) = D50(3,1);
   up_D90(jj,1) = D90(3,1);
   up_eta(jj,1) = eta(3,1);
   up_ent_f(jj,1) = Ei(3,1);
   up_ent_c(jj,1) = Ei(3,7);

   %downstream time variation (M=15)
   dw_eta(jj,1) = eta(15,1);
   dw_D50(jj,1) = D50(15,1);
   dw_D90(jj,1) = D90(15,1);
    
end %for j

%% Plots
% Display results
fprintf("\n total runtime : %g sec \n", Niterations*dt)
fprintf("\n mean elevation : %g m \n", mean(eta))
fprintf("\n max elevation : %g m\n", max(eta))
fprintf("\n min elevation : %g m\n", min(eta))
fprintf("\n slope of levee : %g \n", 8/(sum(eta(:))*dx))
%%
% Levee elevation
figure(3)
    nn = 1:20;
    nn = nn*10;
    plot(nn,eta);  
    title('Sediment surface elevation'); xlabel('levee width [m]'); ylabel('elevation [m]'); axis([0 200 0 2]);
    hold on; %time varying eta
% Aggradation rate at N = 3, 15
figure(4);
    subplot(1,2,1);
    plot(jprint*dt/60,up_eta); axis([0 11000 0 1.8]);
    title('Upstream'); xlabel('T [min]'); ylabel('elevation [m]'); hold on;
    subplot(1,2,2);
    plot(jprint*dt/60,dw_eta); axis([0 11000 0 0.3]);
    title('Downstream'); xlabel('T [min]'); ylabel('elevation [m]'); hold on;
%%
l = length(up_eta);
fprintf("\n aggradation rates : %g, %g \n", (up_eta(l,1)*60)/(jprint(l,1)*dt),(dw_eta(l,1)*60)/(jprint(l,1)*dt))
%%
% Calculating D50 in total deposition
    T_D50 = zeros(M,1);
    xtx = 0.5; % grain size such that xtx % of the material
    u = 0;
    for i=1:M
        for y = 2:1:7
                if T_f(i,y) >= xtx && T_f(i,y-1) < xtx
                  u = y; 
                end
        end
        psix = psi(u-1) + ((psi(u)-psi(u-1))/(T_f(i,u)-T_f(i,u-1)))*(xtx - T_f(i,u-1));  % psi values of x% of material
        T_D50(i,1) = 2^(psix); %grain size values of x% of material [mm]scale
    end
    T_Psibar(i) = sum(T_Pi(i,:)*psi(:)); %mean psi values
    Dm = 2.^(T_Psibar)/1000; %mean grain size (mm)
 
% Calculating D90 in total deposition
    T_D90 = zeros(M,1);
    xtx = 0.9; % grain size such that xtx % of the material
    u = 0;
    for i=1:M
        for y = 2:1:7
                if T_f(i,y) >= xtx && T_f(i,y-1) < xtx
                  u = y; 
                end
        end
        psix = psi(u-1) + ((psi(u)-psi(u-1))/(T_f(i,u)-T_f(i,u-1)))*(xtx - T_f(i,u-1));  % psi values of x% of material
        T_D90(i,1) = 2^(psix); %grain size values of x% of material
    end
%%
% Grainsize Curves across the levee width  
figure(5)
   plot(nn,T_D50); hold on;
   plot(nn,T_D90);
   title('D50 & D90 values at each position'); xlabel('levee width [m]'); ylabel('D [mm]'); axis([0 200 0 0.2]);
   
% Temporal Grainsize variation at N = 3, 15       
figure(6)
   subplot(1,2,1)
   plot(jprint*dt/60,up_D50*1000); hold on; 
   plot(jprint*dt/60,up_D90*1000);
   title('Proximal'); xlabel('time [min]'); ylabel('D [mm]'); axis([0 11000 0.02 0.18]);
   subplot(1,2,2)
   plot(jprint*dt/60,dw_D50*1000); hold on; 
   plot(jprint*dt/60,dw_D90*1000);
   title('Distal'); xlabel('time [min]'); ylabel('D [mm]'); axis([0 11000 0.02 0.18]);
   