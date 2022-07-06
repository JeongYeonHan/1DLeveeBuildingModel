% Plotting of Characteristic Levee Slope and Avulsion Frequency
% With ranges of Overflow Velocity (0.1 ~ 0.5 m/s)
% Median grainsize of SSC (0.05 - 0.25 mm)
% JeongYeon Han and Wonsuck Kim
% Sep 2020
clear all; close all;
%% Assign matrix of Parameters
M=50; %grain size range of inchannel
N=401; %overflow velocity range
tot = M*N; %total matrix
D50 = zeros([tot 1]);
U = zeros([tot 1]);
f = zeros([tot 1]);
S = zeros([tot 1]);
k=0;

%repeatedly run the model with assigned ranges
for i=1:M
    Ps = (i/10)-7.6; %grainsize
    for j=1:N
        Uf=0.1 + 0.001*(j-1); %overflow velocity
        [output1,output2,output3,output4] = AVF(Uf,Ps);
        if output3 > 0 %maintain wedge shape
            continue;
        else
        k=k+1;
        D50(k) = output1; %median grainsize of SSC into floodplain    
        U(k) = Uf; %overflow velocity
        f(k) = 60*60*24*30*12*0.002/output2; %avulsion frequency I_f = 0.002
        S(k) = output4; %characteristic slope
        end
    end
end

%% Plot
%make grid data of slope and avulsion frequency
[xq, yq] = meshgrid(0.5*10^(-4):3*10^(-7):2.501*10^(-4), 0.1:0.001:0.5);
zq = griddata(D50,U,f,xq,yq,'natural');
zq2 = griddata(D50,U,S,xq,yq,'natural');
%% Figure_7A
% Avulsion Frequency Plot
figure(1)
[c,h]=contourf(xq,yq,zq,10000,'EdgeColor','none');
shading interp
c = multigradient('preset','div.cb.BuRd.7');
colormap(c);
hold on;
[c,h]=contour(xq,yq,zq,'k','LevelList', [0.5,1 1.5,2]);
clabel(c,h,[0.5,1 1.5,2],'FontSize',15,'labelspacing',1000);
% Color bar legend
cb1 = colorbar();
set(cb1,'FontSize',15)
title(cb1,'f_A','FontSize',15)
xlabel('d_{50} [m]','FontSize',15)
ylabel('Overflow velocity [m/s]','FontSize',15)
xticks([0.5*10^(-4) 1*10^(-4) 1.5*10^(-4) 2*10^(-4) 2.5*10^(-4)]);
yticks([0.1 0.2 0.3 0.4 0.5]);
% Export plot
print(gcf,'-depsc','-painters','contoursf.eps');

%% Figure_7B
% Characteristic levee slope
figure(2)
[c2,h2]=contourf(xq,yq,zq2,1000,'EdgeColor','none');
shading interp
c2 = multigradient('preset','div.cb.BuRd.7');
colormap(c2);
hold on;
[c2,h2]=contour(xq,yq,zq2,'k','LevelList', [0.1,0.15,0.2,0.25,0.3,0.35,0.4]);
clabel(c2,h2,[0.1,0.15,0.2,0.25,0.3,0.35,0.4],'FontSize',13,'labelspacing',1000);
% Color bar legend
cb1 = colorbar();
set(cb1,'FontSize',15)
title(cb1,'S_r','FontSize',15)
cb1.Ticks = [0.1,0.15,0.2,0.25,0.3,0.35,0.4];
xlabel('d_{50} [m]','FontSize',15)
ylabel('Overflow velocity [m/s]','FontSize',15)
xticks([0.5*10^(-4) 1*10^(-4) 1.5*10^(-4) 2*10^(-4) 2.5*10^(-4)]);
yticks([0.1 0.2 0.3 0.4 0.5]);
% Export plot
print(gcf,'-depsc','-painters','contoursS.eps');

%% 
function [output1, output2, output3, output4] = AVF(UI,Ps)
%% 
g = 9.81; %acceleration of gravity
nu = 1*10^(-6); %kinematic viscosity
ro = 1000; %water density
ros = 2650; %sediment density
R = (ros/ro) -1; %submerged specific gravity
%% model parameters
M = 20; %number of nodes
L = 200; %levee length [m]
%% 
%overflow properties
HI = 4; %initial overflow depth [m]
Ht = 8;
%%
dx = L/M; %step length
dt = 10; %time variation [s]
Niterations = 5000000; %number of total iterations 
%% model matrix
%assign matrix of parameters
qs = zeros(M,7); %sediment discharge
qst = zeros(M,1); %total sediment discharge of each node
cibar = zeros(M,7); % sediment concentration
qsoi = zeros(7,1); %initial sediemnt discharge of each sediment grainsize range
Pi = zeros(M,7); %fractional distribution of ith grain size at each node in the flood flow
Ds = zeros(M,7); %deposit rate for each grain size on one iteration
Dst = zeros (M,1); %total rate at each node on one iteration
%% initial conditions
%assign initial slope of floodplain
eta = zeros (M,1); %bed elevation
SI = 0.0001; %initial slope
etaI = L*SI; %inital elevation at levee front
%setup initial bed surface
for i=1:M
    eta(i) = etaI - SI*dx*i; %initial elevation at each node
end
%% input grain size
psi = Ps:1:(Ps+6);
D = 2.^(psi); %sediment grain size in mm
D = D/1000; %sediment grain size in m
psibar = mean(psi);     %-3; %mean psi
sigma = 0.8;    %standard deviation for normal distribution
p =(1/(sqrt(2*pi)*sigma))*exp(-(((psi-psibar)./sigma).^2)/2); %gaussian grain size distribution
%settling velocity of each grain size using Ferguson and Church (2004)
%the constants, C1 = 18 and C2 = 1 for natural grains
ws=(R*g.*(D.^2))./(18*nu+(0.75*R*g.*D.^3).^0.5);
Repi = ((R*g*D).^(0.5)).*D/nu;  % particle Reynolds Number
qsot = 0.0003;
qsoi = qsot.*p;
f_qso = cumsum(p);
%% initial condition_2
%initial concentration of suspended sediment
alpha = 0; % constant for initial suspended sediment concentration, 0 for fresh water
for i=1:M
    cibar(i,:) = (alpha.*qsoi)./(UI*HI);
end

%calculating D50
xtx = 0.5; % grain size such that xtx % of the material
u = 0;
     for y = 2:1:7
         if f_qso(y) >= xtx && f_qso(y-1) < xtx
            u = y; 
         end
     end
psi_50 = psi(u-1) + ((psi(u)-psi(u-1))/(f_qso(u)-f_qso(u-1)))*(xtx - f_qso(u-1));  % psi values of x% of material
D50_qso = 2^(psi_50)/1000; %grain size values of x% of material
output1 = D50_qso;

%%
    for j=1:1:Niterations
        for i=1:1:M
            for k=1:1:7 %7 grain size
                if i==1
                    cibar(i,k) = cibar(i,k)+ (((qsoi(k)-qs(i,k))/dx)-(ws(k)*cibar(i,k)))*dt/HI; %sediment concentration at first node using the ghost node
                else
                    cibar(i,k) = cibar(i,k)+ (((qs(i-1,k)-qs(i,k))/dx)-(ws(k)*cibar(i,k)))*dt/HI;  %sediment concentration at other nodes
                end
                if cibar(i,k) < 0
                    cibar(i,k) = 0;
                end
                qs(i,k) = cibar(i,k)*UI*HI; %sediment discharge at each node
                Ds(i,k) = ws(k)*cibar(i,k); %deposit rate  for each grain size at each node on one iteration 
            end %for k
                qst(i) = sum(qs(i,:)); %total sediment discharge at each node 
                Dst(i) = sum(Ds(i,:)); %total deposit rate at each node on one iteration
                eta(i) = eta(i) + Dst(i)*dt; %update sediment surface elevation (Exner equation)
                Pi(i,:) = qs(i,:)./qst(i);  % fraction function of each suspended grain size at each node
             
                % using the total depositional rate compared to individual
                % depositional rate to calculate fraction of ith grainsize
                % range in the deposit  
      %print every 100 iteration
        end 
        if eta(1) >= (Ht-HI)
%             fprintf("\n avulsion time : %g sec \n", j*dt);
            break;
        end
    end %for j

%%
output2 = j*dt; %flood occurs for 2 months per year [yr]
output3 = max(diff(eta));
A = sum(eta(:))*dx;
output4 = 8/A; %Hc = 4m
end

   
   