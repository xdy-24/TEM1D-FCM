 function [rho_inv, m_inv, misfit_history, misfit_Bt, centers, U, center_dist_history] = ...
    occam1D_inversion_with_clustering(hpumct, TxCord, recPos, curr, t, obs_data, calType, params)    
%% ----------------------------反演子程序--------------------------------------------------
    %参数变换函数
    function rho = rho_transform(x)
        rho = exp(x);  % 电阻率自然对数变换
    end
    
    function m = m_transform(x)
        m = 1./(1 + exp(-x));  % 充电率logistic变换
    end
    
    function x = rho_inverse_transform(rho)
        x = log(rho);
    end
    
    function x = m_inverse_transform(m)
        x = -log(1./m - 1);
    end

    % 初始化
    n_layers = size(hpumct, 1);
    n_times = length(t);
    
    % 初始化变换后的参数
    ln_rho = rho_inverse_transform(hpumct(:,2));
    m_param = m_inverse_transform(hpumct(:,4));
    
    % 组合优化参数
    opt_params = [ln_rho; m_param];
    n_opt_params = length(opt_params);
    
    % 多起点策略参数
    n_starts =1;
    best_misfit_global = inf;
    best_params_global = opt_params;


    % 参数边界
    rho_bounds = [log(1), log(1000)];  % 在变换空间中的电阻率边界
    m_bounds = [-10, 10];  % 在变换空间中的充电率边界
    
    % 优化参数
    lambda = params.lambda;
    alpha = 1.0;  % 初始步长
    beta = 0.5;   % 步长衰减因子
    
    
%% ------- -聚类参数设置(依据先验信息手动填写)--------
      n_clusters = 3;             % 聚类数量
      target_rho1=rho_inverse_transform(100);%目标电阻率
      target_m1=m_inverse_transform(0.1);%目标充电率
      target_rho2=rho_inverse_transform(200);%目标电阻率2
      target_m2=m_inverse_transform(0.6);%目标充电率2
      back_rho1=rho_inverse_transform(800);%背景电阻率1
      back_m1=m_inverse_transform(0.1);%背景充电率1
      %back_rho2=rho_inverse_transform(800);%背景电阻率2
      %back_m2=m_inverse_transform(0.1);%背景充电率2
      target_centers = [target_rho1, target_m1;target_rho2,target_m2; back_rho1, back_m1];  % 目标中心 (电阻率, 充电率)

%%    
    % 收敛控制参数
    misfit_history = zeros(params.max_iter * n_starts, 1);
    min_improvement = 1e-6;
    stable_window = 3;
    misfit_window = zeros(stable_window, 1);
    
    % 聚类参数初始化
    params_for_clustering = [log(hpumct(:,2)), -log(1./hpumct(:,4) - 1)];
    [centers, U] = fcm(params_for_clustering, n_clusters);
    center_dist_history = zeros(params.max_iter, size(centers, 1));
    
    % 打开日志文件
    log_fid = fopen('inversion_log.txt', 'w');
    bt_log_fid = fopen('bt_data_log.txt', 'w');
    param_log_fid = fopen('parameter_log.txt', 'w');
    
    fprintf(log_fid, '反演迭代日志\n----------------------------------------\n');
    fprintf(log_fid, '迭代次数\t误差\t阻尼因子\t步长\n----------------------------------------\n');
    fprintf(bt_log_fid, '迭代次数\t时间(s)\t预测数据\t观测数据\n----------------------------------------\n');
    fprintf(param_log_fid, 'Iteration\tLayer\tRho\tM\n----------------------------------------\n');

    %残差
    log_filename = 'residuals_log.txt';
    fid = fopen(log_filename, 'w');
    if fid == -1
        error('无法创建日志文件');
    end
    fprintf(fid, '迭代次数\t数据拟合残差\t光滑约束残差\t聚类约束残差\t目标引导残差\n');
    %
    
    % 主迭代循环
    for start = 1:n_starts
        % 随机扰动初始值
        if start > 1
            rand_factor = 0.05;  % 20%的随机扰动
            opt_params = best_params_global .* (1 + rand_factor * (rand(size(opt_params)) - 0.5));
        end
        
        % 存储上一次迭代的参数
        prev_params = opt_params;
        
        % 主迭代循环
        for iter = 1:params.max_iter

         % 权重参数
         progress = iter / params.max_iter;
         w_d = 1.0;   % 数据拟合权重
         w_m = 0.05;  % 平滑约束权重
         w_c = 0.2 * (1 - exp(-3 * progress));  % S型增长,聚类约束
         w_t = 0.2 * (1 - exp(-5 * progress));  % S型增长，目标引导

            % 更新模型参数
            hpumct(:,2) = rho_transform(opt_params(1:n_layers));
            hpumct(:,4) = m_transform(opt_params(n_layers+1:end));
            
            % 记录参数
            for i = 1:n_layers
                fprintf(param_log_fid, '%d\t%d\t%.2f\t%.2f\n', iter, i, hpumct(i,2), hpumct(i,4));
            end
            
            % 正演计算
            pred_data = zeros(n_times, 1);
            for j = 1:n_times
                bt = mexBtFwdArbitraryLoop(hpumct, TxCord, recPos, curr, t(j), calType);
                if ~isempty(bt) && length(bt) >= 3
                    pred_data(j) = bt(3);
                    fprintf(bt_log_fid, '%d\t%.6e\t%.6e\t%.6e\n', iter, t(j), pred_data(j), obs_data(j));
                end
            end
            
            % 误差计算
            [misfit, residuals, misfit_Bt] = calculate_misfit_with_clustering(obs_data, pred_data, hpumct, ...
                target_centers, w_d, w_m, w_c, w_t, n_clusters);
            misfit_history((start-1)*params.max_iter + iter) = misfit;
            misfit_window = [misfit_window(2:end); misfit];
            
            % 聚类分析
            ln_rho = log(hpumct(:,2));  % 归一化后的电阻率
            m_param = -log(1./hpumct(:,4) - 1);  % 归一化后的充电率
            params_for_clustering = [ln_rho, m_param];
            [centers, U] = guided_fcm(params_for_clustering, target_centers, n_clusters);
            
            % 记录聚类中心距离
            for k = 1:n_clusters
                [~, cluster_idx] = max(U(:,k));
                center_dist_history(iter, k) = norm(centers(cluster_idx,:) - target_centers(cluster_idx,:));
            end
            
            % 输出诊断信息
            fprintf('[Iter %3d] Misfit=%.3e | ', iter, misfit);
            fprintf('Cluster Dist: %.2f, %.2f %.2f， %.2f， %.2f| ', center_dist_history(iter,:));
            fprintf('Weights: w_d=%.2f, w_c=%.2f, w_t=%.2f | ', w_d, w_c, w_t);
            fprintf('Lambda=%.2e\n', lambda);
            fprintf(log_fid, '%d\t%.3e\t%.3e\t%.3e\n', iter, misfit, lambda, alpha);
            
            % 动态调整阻尼因子
            if iter > 1
                if misfit < misfit_history((start-1)*params.max_iter + iter-1)
                    lambda = max(lambda/params.v, 1e-6);
                else
                    lambda = min(lambda * 2, 10);
                end
            end
            
            % 自适应步长调整
            if iter > 1
                misfit_ratio = misfit / misfit_history((start-1)*params.max_iter + iter-1);
                if misfit_ratio < 0.8
                    alpha = min(alpha * 1.5, 2.0);
                elseif misfit_ratio > 1.2
                    alpha = max(alpha * 0.5, 0.2);
                end
            end
            
            % 收敛性检查
            if iter >= stable_window
                rel_change = abs(misfit - misfit_window(1)) / misfit_window(1);
                param_change = norm(opt_params - prev_params) / norm(opt_params);
                
               % if rel_change < 1e-4 && param_change < 1e-4 && lambda < 1e-3
               if misfit < 5e-2 || iter >= params.max_iter
                    fprintf('Start %d: 在迭代%d次后收敛\n', start, iter);
                    break;
                end
            end
            
            % 计算雅可比矩阵
            J = calculate_jacobian_with_clustering(hpumct, TxCord, recPos, curr, t, calType, opt_params, ...
                params.d, target_centers, w_d, w_m, w_c, w_t, n_layers, n_clusters);
            
            % 构建正规方程
            H = J' * J;
            g = -J' * residuals;
            
            % 添加阻尼项
            temp = lambda * max(abs(diag(H)));
            H = H + temp * eye(size(H));
            
            % 使用LSQR求解正规方程
            [delta, ~] = lsqr(H, g, 1e-6, 100);
            
            % 步长控制
            max_update = 0.5;
            scale = min(1, max_update / max(abs(delta)));
            delta = scale * delta;
            
            % 线搜索
            alpha_current = alpha;
            best_misfit = inf;
            best_params = opt_params;
            
            for trial = 1:10
                test_params = opt_params + alpha_current * delta;
                
                % 应用参数边界
                test_params(1:n_layers) = max(min(test_params(1:n_layers), rho_bounds(2)), rho_bounds(1));
                test_params(n_layers+1:end) = max(min(test_params(n_layers+1:end), m_bounds(2)), m_bounds(1));
                
                % 更新测试模型
                test_hpumct = hpumct;
                test_hpumct(:,2) = rho_transform(test_params(1:n_layers));
                test_hpumct(:,4) = m_transform(test_params(n_layers+1:end));
                
                % 计算测试模型响应
                test_pred = zeros(n_times,1);
                for j = 1:n_times
                    bt = mexBtFwdArbitraryLoop(test_hpumct, TxCord, recPos, curr, t(j), calType);
                    if ~isempty(bt) && length(bt) >= 3
                        test_pred(j) = bt(3);
                    end
                end
                
                % 计算测试误差
                [test_misfit, ~, ~] = calculate_misfit_with_clustering(obs_data, test_pred, test_hpumct, ...
                    target_centers, w_d, w_m, w_c, w_t, n_clusters);
                
                if test_misfit < best_misfit
                    best_misfit = test_misfit;
                    best_params = test_params;
                end
                
                if test_misfit < 0.98*misfit
                    opt_params = test_params;
                    misfit = test_misfit;
                    alpha = min(alpha_current*1.5, 2.0);
                    break;
                else
                    alpha_current = alpha_current * beta;
                end
            end
            
            % 更新参数
            opt_params = best_params;
            prev_params = opt_params;
            
            % 更新全局最优解
            if misfit < best_misfit_global
                best_misfit_global = misfit;
                best_params_global = opt_params;
            end
             % 计算并记录残差
            data_misfit = calculate_data_misfit(obs_data, pred_data);
            smooth_misfit = calculate_smooth_misfit(hpumct);
            cluster_misfit = calculate_cluster_misfit(hpumct, centers, U);
            target_misfit = calculate_target_misfit(centers, target_centers);
            
            % 记录残差到日志文件
           
            fprintf(fid, '%d\t%.6e\t%.6e\t%.6e\t%.6e\n', iter, data_misfit, smooth_misfit, cluster_misfit, target_misfit);
        end
    end
    
    % 使用全局最优解
    opt_params = best_params_global;

   
            fclose(fid);

    % 关闭日志文件
    fclose(log_fid);
    fclose(bt_log_fid);
    fclose(param_log_fid);
    
    % 最终结果
    rho_inv = rho_transform(opt_params(1:n_layers));
    m_inv = m_transform(opt_params(n_layers+1:end));
    misfit_history = misfit_history(1:iter);

    % 分类层
    [is_polarized, ~] = classify_layers(hpumct, centers, U, target_centers);
 
   %% -----------------------------------------图形绘制-------------------------------------- 
% 绘制聚类结果
    figure;
    scatter(log(rho_inv), -log(1./m_inv - 1), 100, U(2,:), 'filled');
    hold on;
    plot(target_centers(:,1), target_centers(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
    xlabel('Resistivity','FontSize', 15);
    ylabel('m','FontSize', 15);
    title('FCM','FontSize', 15);
    colorbar;
    %grid on;
    


% 绘制观测数据与预测数据对比
figure;
loglog(t, abs(obs_data), 'k*', 'MarkerSize', 8); hold on;
loglog(t, abs(pred_data), 'r-', 'LineWidth', 1.5);
legend('obs', 'pre','FontSize', 13);
xlabel('t（s）','FontSize', 15);
ylabel('dBz/dt(T/s)','FontSize', 15);
%title('观测与预测数据对比');
%grid on;

% 绘制相对误差图
figure;
valid_idx = obs_data ~= 0;
relative_error = abs((pred_data(valid_idx) - obs_data(valid_idx)) ./ obs_data(valid_idx));
semilogx(t(valid_idx), relative_error*100, 'b-o', 'LineWidth', 1.5);
xlabel('t（s）','FontSize', 15);
ylabel('Relative error（%）','FontSize', 15);
%title('时间道相对误差');
%grid on;

% 读取日志文件
log_data = readtable(log_filename);

% 绘制残差曲线
% figure;

% % 数据拟合残差
% subplot(2,2,1);
% plot(log_data{:, 1}, log_data{:, 2}, 'b-');
% title('数据拟合残差');
% xlabel('迭代次数');
% ylabel('残差');
% %grid on;
% 
% % 光滑约束残差
% subplot(2,2,2);
% plot(log_data{:, 1}, log_data{:, 3}, 'g-');
% title('光滑约束残差');
% xlabel('迭代次数');
% ylabel('残差');
% grid on;
% 
% % 聚类约束残差
% subplot(2,2,3);
% plot(log_data{:, 1}, log_data{:, 4}, 'r-');
% title('聚类约束残差');
% xlabel('迭代次数');
% ylabel('残差');
% %grid on;
% 
% % 目标引导残差
% subplot(2,2,4);
% plot(log_data{:, 1}, log_data{:, 5}, 'k-');
% title('目标引导残差');
% xlabel('迭代次数');
% ylabel('残差');
% %grid on;

end

function misfit = calculate_data_misfit(obs_data, pred_data)
    misfit = norm(abs((obs_data - pred_data)/obs_data));
end

function misfit = calculate_smooth_misfit(hpumct)
    n_layers = size(hpumct, 1);
    misfit = 0;
    for i = 1:n_layers-1
        misfit = misfit + (log(hpumct(i+1,2)) - log(hpumct(i,2)))^2;
        misfit = misfit + (hpumct(i+1,4) - hpumct(i,4))^2;
    end
end

function misfit = calculate_cluster_misfit(hpumct, centers, U)
    n_layers = size(hpumct, 1);
    n_clusters = size(centers, 1);
    misfit = 0;
    for i = 1:n_layers
        for k = 1:n_clusters
            misfit = misfit + U(k,i)^2 * norm([log(hpumct(i,2)), hpumct(i,4)] - centers(k,:))^2;
        end
    end
end

function misfit = calculate_target_misfit(centers, target_centers)
    n_clusters = size(centers, 1);
    misfit = 0;
    for k = 1:n_clusters
        misfit = misfit + norm(centers(k,:) - target_centers(k,:))^2;
    end
end