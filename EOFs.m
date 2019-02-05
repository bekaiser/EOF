% 12.805 homework 6 (EOFs)
% Bryan Kaiser
% 5/4/16

close all; clear all; clc;
addpath('./functions'); run('./functions/mcolormaps.m');

% load data
data  = importdata('code2c3lp.mat');
t = data(:,1); % days, from 1992
u = data(:,2:6)'; % cm/s, (onshore) ocean current velocity
v = data(:,7:11)'; % cm/s, (alongshore) ocean current velocity
uw = data(:,12); % m/s, (onshore) wind velocity
vw = data(:,13); % m/s, (alongshore) wind velocity
h = -[5 10 15 35 70]; % m, depths
[H,T] = meshgrid(h,t);

[M,N] = size(u); % M = spatial resolution, N = temporal resolution

% remove the time mean at each location:
u = u-mean(mean(u)); v = v-mean(mean(v)); 

% covariance & singular value decomposition:
C_uu = cov(u',1); [U1,S1,V1T] = svd(u); 
C_vv = cov(v',1); [U2,S2,V2T] = svd(v); 

% EOF by covariance method:
EOF_cov_u = C_uu*U1;
EOF_cov_v = C_vv*U2;

% EOF by svd method:
gamma_u = (S1*S1')/N; % eigenvalues on the diagonal 
gamma_v = (S2*S2')/N;
EOF_svd_u = U1*gamma_u;
EOF_svd_v = U2*gamma_v;

% check that the methods are the same:
sum(sum(EOF_cov_u-EOF_svd_u)); % ~ machine precision, checks
sum(sum(EOF_cov_v-EOF_svd_v)); % ~ machine precision, checks

% first spatial EOF, structure:
EOF_1u = U1(:,1);
EOF_1v = U2(:,1);
figure;
p1 = plot(EOF_1u,h,'Color',[0.2 0.2 0.8]); hold on
p2 = plot(EOF_1v,h,'Color',[0.2 0.8 0.2]); 
xlabel('EOF1'); ylabel('z (m)'); hl=legend([p1 p2],'u','v'); 
set(gca,'FontSize',12); set(hl,'FontSize',14);

% variance:
eig_u = diag(gamma_u); % eigenvalues on the diagonal of S
eig_v = diag(gamma_v); 
var_u = eig_u./sum(eig_u); % total variance = sum of eigenvalues
var_v = eig_v./sum(eig_v);
figure;
p1 = plot(1:length(eig_u),var_u,'Color',[0.2 0.2 0.8]); hold on
p2 = plot(1:length(eig_v),var_v,'Color',[0.2 0.8 0.2]); 
xlabel('EOF','FontSize',14); ylabel('% variance','FontSize',14);
hl=legend([p1 p2],'u','v'); set(gca,'FontSize',12); set(hl,'FontSize',14);

% hovmuller plot, u
figure; surf(T,H,u'); view(2); %shading flat 
colorbar; colormap(redblue/255); %colormap(flipud(colormap));  
shading flat; grid off; axis tight
xlabel('t (days)','FontSize',14); ylabel('z (m)','FontSize',14); 
title('u (cm/s)','FontSize',14); set(gca,'FontSize',14); caxis([-25,25]);

% hovmuller plot,
figure; surf(T,H,v'); view(2); %shading flat 
colorbar; colormap(redblue/255); %colormap(flipud(colormap));  
shading flat; grid off; axis tight
xlabel('t (days)','FontSize',14); ylabel('z (m)','FontSize',14); 
title('v (cm/s)','FontSize',14); set(gca,'FontSize',14); caxis([-55,55]);

% expansion coefficients:
b1 = U1'*u; b2 = U2'*v;

% u expansion coefficient plots
figure;
p1 = plot([0:length(b1(1,:))-1],b1(1,:),'Color',[0.2 0.2 0.8]); hold on
p2 = plot([0:length(b1(1,:))-1],b1(2,:),'Color',[0.2 0.4 0.6]); 
p3 = plot([0:length(b1(1,:))-1],b1(3,:),'Color',[0.2 0.5 0.5]); 
p4 = plot([0:length(b1(1,:))-1],b1(4,:),'Color',[0.2 0.6 0.4]); 
p5 = plot([0:length(b1(1,:))-1],b1(5,:),'Color',[0.2 0.8 0.2]); 
xlabel('days','FontSize',14); ylabel('b_i(t)','FontSize',14); title('u','FontSize',14);
hl = legend([p1 p2 p3 p4 p5],'1','2','3','4','5'); set(gca,'FontSize',14)
axis([0,length(b1(1,:))-1,-50,60]); set(hl,'FontSize',14);

% v expansion coefficient plots
figure;
p1 = plot([0:length(b2(1,:))-1],b2(1,:),'Color',[0.2 0.2 0.8]); hold on
p2 = plot([0:length(b2(1,:))-1],b2(2,:),'Color',[0.2 0.4 0.6]); 
p3 = plot([0:length(b2(1,:))-1],b2(3,:),'Color',[0.2 0.5 0.5]); 
p4 = plot([0:length(b2(1,:))-1],b2(4,:),'Color',[0.2 0.6 0.4]); 
p5 = plot([0:length(b2(1,:))-1],b2(5,:),'Color',[0.2 0.8 0.2]); 
xlabel('days','FontSize',14); ylabel('b_i(t)','FontSize',14); title('v','FontSize',14)
hl=legend([p1 p2 p3 p4 p5],'1','2','3','4','5'); set(gca,'FontSize',14);
axis([0,length(b1(1,:))-1,-100,150]); set(hl,'FontSize',14);

% reconstruction the series without EOFs 3-5:
U1(:,3:5) = zeros(5,3); U2(:,3:5) = zeros(5,3); 
ufilt = U1*b1; vfilt = U2*b2;

figure;
p1 = plot(t,v(1,:),'b','LineWidth',1.2); hold on
p2 = plot(t,vfilt(1,:),'r','LineWidth',1.2); 
p3 = plot(t,vfilt(1,:)-v(1,:),'k','LineWidth',1.2); 
xlabel('days','FontSize',14); ylabel('cm/s','FontSize',14);
hl = legend([p1 p2 p3],'v','v_{1,2}','v-v_{1,2}'); set(hl,'FontSize',14);
set(gca,'FontSize',12);
axis([t(1),t(end),-70,70])


%==============================================================================
% problem 2: objective mapping
% use objective mapping to estimate Newport sea level for the first 25 days
% of data for dt = 0.1 days. 

% load data
data = importdata('hwa4ar.mat');
y = data(1:25,3); % cm, Newport sea level (size [1 x 25])
t1 = data(1:25,1); % days, dt = 1 day
t=[1464:0.1:1489]; % days (dt = 0.1 day intervals, size [1 x 251])
M = length(y); % number of observations = 25
N = length(t); % number of solutions sought = 251

% 1) covariance matrix as a function of lag
tau = [0:0.1:25];
Rx = exp(-(tau./4).^2).*cos(tau.*(2*pi/8)).*100;
Rxx = toeplitz(Rx); % size [251 x 251]
% checks:
%sum(Rxx(1,:)-Rx)
%size(Rxx)
%sum(diag(Rxx)./100)

% 2) second moment of the observational error (size [25 x 25])
Rnn = eye(M).*10;

% 3) make the observational mapping matrix E (must be size [25 x 251])
E = zeros(M,N); locs = [6:10:251];
for m = 1:M
    E(m,locs(m)) = 1;
end

% 4) Rxy and Ryy
Rxy = Rxx*E'; Ryy = E*Rxx*E'+Rnn;

% 5) Objectively mapped estimate
x = Rxy/Ryy*y;

% 6) Uncertainty
P = Rxx-Rxx*E'/Ryy*E*Rxx; Perr = (diag(P)).^(1/2);

figure; 
p2 = plot(t,x,'Color',[0.2 0.2 0.8],'LineWidth',1.2); hold on 
p3 = plot(t,x+Perr,':','Color',[0.2 0.2 0.8],'LineWidth',1.2);
p3 = plot(t,x-Perr,':','Color',[0.2 0.2 0.8],'LineWidth',1.2);
p1 = plot(t1,y,'o','Color',[0.3 0.7 0.3],'LineWidth',1.2); 
xlabel('days','FontSize',14); ylabel('ssh (cm)','FontSize',14); 
hl=legend([p1 p2 p3],'y','x','error'); set(hl,'FontSize',14); 
axis([t(1),t(end),-25,20]); set(gca,'FontSize',12)




