function J = calculate_jacobian_with_clustering(hpumct, TxCord, recPos, curr, t, calType, opt_params, d, ...
    target_centers, w_d, w_m, w_c, w_t, n_layers, n_clusters)
    
    % 初始化
    t = t(:);
    n_times = length(t);
    n_params = length(opt_params);
    
    % 计算约束数量
    n_model_constraints = n_layers - 1;
    n_cluster_constraints = n_layers;
    total_constraints = n_times + n_model_constraints + n_cluster_constraints;
    
    % 初始化Jacobian矩阵
    J = zeros(total_constraints, n_params);
    
    % 带符号的对数变换函数
    function y = signed_log(x)
        y = sign(x) .* log(1 + abs(x));
    end
    
    % 1. 数据域Jacobian
    for i = 1:n_params
        % 计算扰动量
        if i <= n_layers  % 电阻率参数
            delta = d * abs(opt_params(i));
        else  % 充电率参数
            delta = d * 10.0;  % 增加充电率扰动量
        end
        
        hpumct_plus = hpumct;
        hpumct_minus = hpumct;
        
        if i <= n_layers
            hpumct_plus(i,2) = exp(opt_params(i) + delta);
            hpumct_minus(i,2) = exp(opt_params(i) - delta);
        else
            idx = i - n_layers;
            hpumct_plus(idx,4) = opt_params(i) + delta;
            hpumct_minus(idx,4) = opt_params(i) - delta;
        end
        
        % 计算正演响应差分
        for j = 1:n_times
            try
                % 计算当前模型、正扰动模型和负扰动模型的正演响应
                bt_current = mexBtFwdArbitraryLoop(hpumct, TxCord, recPos, curr, t(j), calType);
                bt_plus = mexBtFwdArbitraryLoop(hpumct_plus, TxCord, recPos, curr, t(j), calType);
                bt_minus = mexBtFwdArbitraryLoop(hpumct_minus, TxCord, recPos, curr, t(j), calType);
                
                if ~isempty(bt_plus) && ~isempty(bt_minus) && ~isempty(bt_current) && ...
                   length(bt_plus) >= 3 && length(bt_minus) >= 3 && length(bt_current) >= 3
                    % 计算相对变化率
                    rel_plus = (bt_plus(3) - bt_current(3)) / abs(bt_current(3));
                    rel_minus = (bt_minus(3) - bt_current(3)) / abs(bt_current(3));
                    
                    % 对相对变化进行带符号的对数变换
                    log_plus = signed_log(rel_plus);
                    log_minus = signed_log(rel_minus);
                    
                    % 计算导数
                    deriv = (log_plus - log_minus) / (2*delta);
                    
                    % 参数敏感度调整
                    if i <= n_layers
                        sensitivity = 1.0;  % 电阻率敏感度
                    else
                        sensitivity = 20.0;  % 增加充电率敏感度
                    end
                    
                    J(j,i) = w_d * deriv * sensitivity;
                end
            catch
                J(j,i) = 0;
            end
        end
    end
    
    % 2. 模型平滑约束的Jacobian
    for i = 1:n_layers-1
        row_idx = n_times + i;
        
        if i <= n_layers-1
            % 电阻率平滑约束
            J(row_idx,i) = -w_m;
            J(row_idx,i+1) = w_m;
        end
        
        % 充电率平滑约束
        if i+n_layers <= n_params
            % 增加充电率约束的权重
            w_m_charge = w_m * 2.0;
            J(row_idx,i+n_layers) = -w_m_charge;
            J(row_idx,i+1+n_layers) = w_m_charge;
        end
    end
    
    % 3. 聚类约束Jacobian
    ln_rho = log(hpumct(:,2));  % 归一化后的电阻率
    m_param = -log(1./hpumct(:,4) - 1);  % 归一化后的充电率
    params_for_clustering = [ln_rho, m_param];
    [centers, U] = guided_fcm(params_for_clustering, target_centers, n_clusters);
    
    for i = 1:n_layers
        row_idx = n_times + n_layers - 1 + i;
        
        for j = 1:n_params
            % 计算扰动后的参数
            hpumct_plus = hpumct;
            if j <= n_layers
                delta = d * abs(opt_params(j));
                hpumct_plus(j,2) = exp(opt_params(j) + delta);
            else
                idx = j - n_layers;
                delta = d;
                hpumct_plus(idx,4) = opt_params(j) + delta;
            end
            
            % 计算原始和扰动后的聚类参数
            log_rho_orig = log(hpumct(:,2));
            log_rho_plus = log(hpumct_plus(:,2));
            
            m_orig = hpumct(:,4);
            m_plus = hpumct_plus(:,4);
            
            params_orig = [log_rho_orig, m_orig];
            params_plus = [log_rho_plus, m_plus];
            
            % 确保数据有效
            params_orig(isnan(params_orig)) = 0;
            params_orig(isinf(params_orig)) = 0;
            params_plus(isnan(params_plus)) = 0;
            params_plus(isinf(params_plus)) = 0;
            
            % 使用Guided FCM计算
            [centers_orig, U_orig] = guided_fcm(params_orig, target_centers, n_clusters);
            [centers_plus, U_plus] = guided_fcm(params_plus, target_centers, n_clusters);
            
            [max_u_orig, cluster_idx_orig] = max(U_orig(:,i));
            [max_u_plus, cluster_idx_plus] = max(U_plus(:,i));
            
            % 计算距离变化
            dist_orig = max_u_orig^2 * norm(params_orig(i,:) - centers_orig(cluster_idx_orig,:));
            dist_plus = max_u_plus^2 * norm(params_plus(i,:) - centers_plus(cluster_idx_plus,:));
            
            target_dist_orig = norm(centers_orig(cluster_idx_orig,:) - target_centers(cluster_idx_orig,:));
            target_dist_plus = norm(centers_plus(cluster_idx_plus,:) - target_centers(cluster_idx_plus,:));
            
            % 计算导数
            J(row_idx,j) = w_c * (dist_plus - dist_orig) / delta + w_t * (target_dist_plus - target_dist_orig) / delta;
        end
    end
    
    % 数值稳定性处理
    J(isnan(J)) = 0;
    J(isinf(J)) = 0;
    
    % 归一化Jacobian矩阵
    norm_factor = max(abs(J(:))) + eps;
    J = J / norm_factor;
end 