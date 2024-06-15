% Author : Cheolmin Jeong
% E-mail : cmin87394@gmail.com
% Reference URL : https://www.dbpia.co.kr/journal/articleDetail?nodeId=NODE10607318


close all; clear all; clc;

%% Electric Power Steering System
% real system
[A, B, C, D] = SystemModel;

Tc = 0.01;
T = 0:Tc:10;
len = length(T);

var = length(C{1});
x = zeros(var,len);
y = zeros(1,len); 

Q = diag([1 0.7 0.5 0.1]);
Qc = chol(Q);
w = [Qc(1,1) Qc(2,2) Qc(3,3) Qc(4,4)]';

R = 1;
v = chol(R);

for k = 1:len-1
    u(k) = sin(0.01*pi*k);
    if k <= 500
        x(:,k+1) = A{1}*x(:,k) + B{1}*u(k) + w*(2*rand-1);
        y(k+1) = C{1}*x(:,k+1) + v*(2*rand-1);
    else
        x(:,k+1) = A{2}*x(:,k) + B{2}*u(k) + w*(2*rand-1);
        y(k+1) = C{2}*x(:,k+1) + v*(2*rand-1);
    end
end
u(len) = sin(0.01*pi*len);

% ideal system(no noise)
x_tr = zeros(var,len);
y_tr = zeros(1,len);
for k = 1:len-1
    if k <= 500
        x_tr(:,k+1) = A{1}*x_tr(:,k) + B{1}*u(k);
        y_tr(k+1) = C{1}*x_tr(:,k+1);
    else
        x_tr(:,k+1) = A{2}*x_tr(:,k) + B{2}*u(k);
        y_tr(k+1) = C{2}*x_tr(:,k+1);
    end
end

%% Interacting Multiple Model
% Number of Models
m = 2;

% Initial Probability
mu_i = [0.5 0.5];
p_ij = [0.99 0.01;
        0.01 0.99];

x_h{1} = zeros(var,1);
P_h{1} = diag([1 1 1 1]);
x_h{2} = zeros(var,1);
P_h{2} = diag([1 1 1 1]);

% Iteration
for k = 1:len
    %% Mixing
    % Normalizing factors for mixing probabilities
    c_j = zeros(1,m);
    for j = 1:m
        for i = 1:m
            c_j(j) = c_j(j) + p_ij(i,j).*mu_i(i);
        end
    end

    % Mixing probabilities
    mu_ij = zeros(m,m);
    for i = 1:m
        for j = 1:m
            mu_ij(i,j) = p_ij(i,j) * mu_i(i) / c_j(j);
        end
    end

    % Calculate the mixed state mean for each filter
    x_0j = cell(1,m);
    for j = 1:m
        x_0j{j} = zeros(var,1);
        for i = 1:m
            x_0j{j} = x_0j{j} + x_h{i}*mu_ij(i,j);
        end
    end

    % Calculate the mixed state covariance for each filter
    P_0j = cell(1,m);
    for j = 1:m
        P_0j{j} = zeros(var,var);
        for i = 1:m
            P_0j{j} = P_0j{j} + mu_ij(i,j)*(P_h{i} + (x_h{i}-x_0j{j})*(x_h{i}-x_0j{j})');
        end
    end

    %% Prediction
    x_bar = cell(1,m);
    P_bar = cell(1,m);
    for j = 1:m
        x_bar{j} = A{j}*x_0j{j} + B{j}*u(k);
        P_bar{j} = A{j}*P_0j{j}*A{j}' + Q;
    end
        
    %% Correction
    LH = zeros(1,m);
    for j = 1:m
        S = C{j}*P_bar{j}*C{j}' + R;
        K = P_bar{j}*C{j}'*inv(S);
        P_h{j} = P_bar{j} - K*C{j}*P_bar{j};
        x_h{j} = x_bar{j} + K*y(k) - K*C{j}*x_bar{j};

        IM = C{j}*x_bar{j};
        LH(j) = Gauss_PDF(y(k),IM,S);
    end
    
    %% Merging
    % Calculate the model probabilities
    mu_i = c_j.*LH/sum(c_j.*LH);
    
    % Output the combined updated state mean and covariance
    % Space for estimates
    x_est = zeros(var,1);
    P_est = zeros(var,var);
    
    % Updated state mean
    for j = 1:m
        x_est = x_est + mu_i(j)*x_h{j};
    end
    % Updated state covariance
    for j = 1:m
        P_est = P_est + mu_i(j)*(P_h{j} + (x_h{j}-x_est)*(x_h{j}-x_est)');
    end
    
    y_est = C{1}*x_est;
    
    MP(:,k) = mu_i;
    x_estT(:,k) = x_est;
    y_estT(:,k) = y_est;
    
end

%% RMSE
rms_x1 = sqrt(sum((x(1,:) - x_est(1,:)).^2)/len)
rms_x2 = sqrt(sum((x(2,:) - x_est(2,:)).^2)/len)
rms_x3 = sqrt(sum((x(3,:) - x_est(3,:)).^2)/len)
rms_x4 = sqrt(sum((x(4,:) - x_est(4,:)).^2)/len)
rms_y = sqrt(sum((y(:) - y_est(:)).^2)/len)

%% Plotting
% States
figure(1)
subplot(411)
plot(T,x(1,:),'color','#A4A3A3','linewidth',2)
hold on
plot(T,x_estT(1,:),'color','#FA902D','linewidth',1)
legend('real system', 'IMM','fontsize', 12,'fontweight','bold')
xlabel('time','fontsize',12,'fontweight','bold')
ylabel('x1','fontsize',12,'fontweight','bold')
subplot(412)
plot(T,x(2,:),'color','#A4A3A3','linewidth',2)
hold on
plot(T,x_estT(2,:),'color','#FA902D','linewidth',1)
legend('real system', 'IMM','fontsize', 12,'fontweight','bold')
xlabel('time','fontsize', 12,'fontweight','bold')
ylabel('x2','fontsize', 12,'fontweight','bold')
subplot(413)
plot(T,x(3,:),'color','#A4A3A3','linewidth',2)
hold on
plot(T,x_estT(3,:),'color','#FA902D','linewidth',1)
legend('real system', 'IMM','fontsize', 12,'fontweight','bold')
xlabel('time','fontsize', 12,'fontweight','bold')
ylabel('x3','fontsize', 12,'fontweight','bold')
subplot(414)
plot(T,x(4,:),'color','#A4A3A3','linewidth',2)
hold on
plot(T,x_estT(4,:),'color','#FA902D','linewidth',1)
legend('real system', 'IMM','fontsize', 12,'fontweight','bold')
xlabel('time','fontsize', 12,'fontweight','bold')
ylabel('x4','fontsize', 12,'fontweight','bold')

% Output
figure(2)
plot(T,y(:),'color','#A4A3A3','linewidth',2)
hold on
plot(T,y_estT(:),'color','#FA902D','linewidth',1)
legend('real system', 'IMM','fontsize', 12,'fontweight','bold')
xlabel('time','fontsize', 12,'fontweight','bold')
ylabel('y','fontsize', 12,'fontweight','bold')

% Input
figure(3)
plot(T,u(:))

% Mode probability
figure(4)
plot(T,MP(1,:),'linewidth',1.5)
hold on
plot(T,MP(2,:),'linewidth',1.5)
legend('Model1', 'Model2','fontsize', 12,'fontweight','bold')

