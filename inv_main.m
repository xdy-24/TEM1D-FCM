%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        程序功能：考虑激发极化特征的瞬变电磁聚类约束一维反演       %
%          作    者：谢东阳 
%          邮    箱：2024126079@chd.edu.cn
%          完成时间：2025年4月16日                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% --------------------------模型参数设计---------------------------------
clc
clear ;
%计时
tic; 
% 设置时间道
nt =100; 
% 初始模型参数
thk = [50 50 50 5000]';     % 层厚度
rho = [150 150 150 150]';    % 初始电阻率
miu = ones(4, 1);           % 磁导率
m = [0.3 0.3 0.3 0.3]';           % 充电率
c = [0 0 0.5 0 ]';           % 频率依赖系数
tau = [0 0 0.01 0 ]';        % 时间常数

% 组合地层参数
hpumct = [thk, rho, miu, m, c, tau];

% 真实地层模型参数
rho_true = [100 200 800 ]';    % 真实电阻率
m_true = [ 0.1 0.6 0.1]';    % 真实充电率
hpumct_true = [100 100; 50 200;5000 800];  % 真实模型层厚度和电阻率

% 发射源和接收点参数
Tx = [200 200];
lxhf = 0.5*Tx(1);
lyhf = 0.5*Tx(2);   
TxCord = [lxhf,-lyhf; lxhf,lyhf; -lxhf,lyhf; -lxhf,-lyhf; lxhf,-lyhf];
recPos = [0 0 0];          % 接收点位置
curr = 10.0;                 % 发射电流
calType = 0;                % 计算类型(1-磁场,0-衰减电压)

%% 读取观测数据
fid1 = fopen('realmodel.txt', 'r');
data = fscanf(fid1, '%f', [2, nt]);
t = data(1, :)';
obs_data = data(2, :)';
fclose(fid1);

%% 反演参数设置
params = struct();
params.max_iter =30;       % 最大迭代次数
params.lambda = 1;         % 初始阻尼因子
params.eps = 1e-6;          % 收敛阈值
params.d = 1e-2;            % 求导步长
params.v = 2;               % lambda调节因子
params.target_misfit = 1e-3;% 目标误差值

%%
% 启动并行池
if isempty(gcp('nocreate'))
    parpool('local', 4); % 启动4个线程
end
%% 运行反演
%固定权重版本
[rho_inv, m_inv, misfit_history, misfit_Bt, centers, U, center_dist_history] = ...
    occam1D_inversion_with_clustering(hpumct, TxCord, recPos, curr, t, obs_data, calType, params);

%下面版本A_occam1D_inversion_with_clustering的是预热机制动态调整权重版本
%[rho_inv, m_inv, misfit_history, misfit_Bt, centers, U, center_dist_history] = ...
%    A_occam1D_inversion_with_clustering(hpumct, TxCord, recPos, curr, t, obs_data, calType, params);

% ========== 【2】迭代结束，停止计时并计算耗时 ==========
total_time = toc;  
fprintf('\n==========================================\n');
fprintf('         反演迭代完成！\n');
fprintf(' 总运行时间：%.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);
fprintf('==========================================\n\n');

%% 把“层厚+层值”变成阶梯剖面
function [z, val] = layer2step(h, val_vec)
% h:        n×1 厚度
% val_vec:  n×1 层值（电阻率或充电率）
% z, val:   2n×1 阶梯节点
    n = numel(h);
    z   = zeros(2*n,1);
    val = zeros(2*n,1);
    z(1) = 0;  val(1) = val_vec(1);
    for k = 1:n
        z(2*k-1)   = sum(h(1:k-1));      % 层顶
        z(2*k)     = sum(h(1:k));        % 层底
        val(2*k-1) = val_vec(k);
        val(2*k)   = val_vec(k);
    end
end

%% ------------------真实模型-----------------------------
h_true   = hpumct_true(:,1);          % 5×1 厚度
[z_true, rho_true_step] = layer2step(h_true, rho_true);
[~,      m_true_step]   = layer2step(h_true, m_true);

%% ------------------反演模型-----------------------------
h_inv   = hpumct(:,1);                % 8×1 厚度
rho_inv = rho_inv(:);                 % 8×1 反演电阻率
m_inv   = m_inv(:);                   % 8×1 反演充电率
[z_inv, rho_inv_step]  = layer2step(h_inv, rho_inv);
[~,      m_inv_step]    = layer2step(h_inv, m_inv);

%% ------------------绘图-----------------------------
figure('Position',[100,100,800,600]);

% 电阻率-深度
subplot(1,2,1);
semilogx(rho_true_step, z_true, 'r-', 'LineWidth',2); hold on;
semilogx(rho_inv_step,  z_inv,  'b--','LineWidth',2);
set(gca,'YDir','reverse');
xlabel('Resistivity (Ω·m)','FontSize',15);
ylabel('Depth (m)','FontSize',15);
%title('电阻率-深度剖面','FontSize',14);
legend('True','inversion','FontSize', 13,'Location','southwest');
ylim([0 500]);%-------------限制了深度可根据需要手动调整
xlim([min([rho_true_step; rho_inv_step])*0.8 ...
      max([rho_true_step; rho_inv_step])*1.2]);

% 充电率-深度
subplot(1,2,2);
plot(m_true_step, z_true, 'r-', 'LineWidth',2); hold on;
plot(m_inv_step,  z_inv,  'b--','LineWidth',2);
set(gca,'YDir','reverse');
xlabel('m','FontSize',15);
ylabel('Depth (m)','FontSize',15);
%title('充电率-深度剖面','FontSize',14);
legend('True','inversion','FontSize', 13,'Location','southwest');
ylim([0 500]);%-------------限制了深度可根据需要手动调整
xlim([0 1]);

% 绘制误差收敛曲线
figure;
plot(misfit_history * 100, 'k-o', 'LineWidth', 1.5);  % 乘以100变成百分比
xlabel('Iteration count','FontSize', 15);
ylabel('error（%）','FontSize', 15);
%title('误差收敛曲线（线性尺度）');
%grid on;

% 保存结果
save('inversion_result.mat', 'rho_inv', 'm_inv', 'misfit_history');
%% 关闭并行池
if isempty(gcp('nocreate')) == 0
    delete(gcp);
    disp('并行池已关闭');
end 