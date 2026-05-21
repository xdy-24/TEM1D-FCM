function [misfit, residuals, misfit_Bt] = calculate_misfit_with_clustering(obs_data, pred_data, hpumct, ...
    target_centers, w_d, w_m, w_c, w_t, n_clusters)
    
    % 初始化
    nt = length(obs_data);
    n_layers = size(hpumct,1);
    
    % 计算残差数量
    n_data_residuals = nt;
    n_model_residuals = n_layers - 1;
    n_cluster_residuals = n_layers;
    total_residuals = n_data_residuals + n_model_residuals + n_cluster_residuals;
    residuals = zeros(total_residuals, 1);
    
    % 带符号的对数变换函数
    function y = signed_log(x)
        y = sign(x) .* log(1 + abs(x));
    end
    
    % 1. 数据拟合残差
    valid_idx = obs_data ~= 0;
    if any(valid_idx)
        % 计算相对误差
        rel_error = (pred_data(valid_idx) - obs_data(valid_idx)) ./ abs(obs_data(valid_idx));
        
        % 对相对误差进行带符号的对数变换
        log_error = signed_log(rel_error);
        
        % 时间权重 - 后期时间道更重要
        error_norm = abs(rel_error) / max(abs(rel_error));
        
        % 调整权重：误差较大的时间道分配更大的权重
        time_weight = 1.5 + error_norm;
        
        % 残差赋值
        residuals(1:sum(valid_idx)) = w_d * log_error .* time_weight;
       
        misfit_Bt = abs(rel_error);
    else
        misfit_Bt = zeros(nt, 1);
    end
    
    % 2. 模型平滑约束
    for i = 1:n_layers-1
        idx = n_data_residuals + i;
        
        % 电阻率梯度（对数空间）
        rho_diff = log(hpumct(i+1,2)/hpumct(i,2));
        
        % 充电率梯度（增加权重）
        m_diff = (hpumct(i+1,4) *1- hpumct(i,4)) * 2;  % 增加充电率约束权重
        
        residuals(idx) = w_m * (rho_diff^2 + m_diff^2);
    end
    
    % 3. 聚类约束
    % 使用归一化后的电阻率和充电率
    ln_rho = log(hpumct(:,2));  % 归一化后的电阻率
    m_param = -log(1./hpumct(:,4) - 1);  % 归一化后的充电率
    params_for_clustering = [ln_rho, m_param];
    [centers, U] = guided_fcm(params_for_clustering, target_centers, n_clusters);
    
    for i = 1:n_layers
        idx = n_data_residuals + n_model_residuals + i;
        [max_u, cluster_idx] = max(U(:,i));
        
        % 计算到聚类中心的距离
        dist_to_center = norm(params_for_clustering(i,:) - centers(cluster_idx,:));
        
        % 计算聚类中心到目标中心的距离
        dist_to_target = norm(centers(cluster_idx,:) - target_centers(cluster_idx,:));
        
        % 组合残差
        residuals(idx) = w_c * max_u^2 * dist_to_center + w_t * dist_to_target;
    end
    
    % 计算总体误差
    valid_residuals = residuals(~isnan(residuals) & ~isinf(residuals));
    if isempty(valid_residuals)
        misfit = 1e10;
    else
        misfit = sqrt(mean(valid_residuals.^2));
    end
end