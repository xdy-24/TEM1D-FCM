 function [centers, U] = guided_fcm(data, target_centers, n_clusters)
    % Guided FCM聚类
    % data: 输入数据
    % target_centers: 目标中心
    % n_clusters: 聚类数量
    
    % 初始化隶属度矩阵
    n_data = size(data, 1);
    U = rand(n_clusters, n_data);
    U = U ./ repmat(sum(U, 1), n_clusters, 1);
    
    max_iter = 100;
    tol = 1e-5;
    
    for iter = 1:max_iter
        % 计算聚类中心
        centers = (data' * U')' ./ sum(U, 2);
        
        % 引导项：将聚类中心拉向目标中心
        centers = centers + 0.4* (target_centers - centers);%%%%先验知识权重应该更大
        
        % 更新隶属度矩阵
        dist = zeros(n_clusters, n_data);
        for k = 1:n_clusters
            dist(k, :) = sum((data - repmat(centers(k,:), n_data, 1)).^2, 2);
        end
        dist = dist + eps;  % 避免除零
        
        U_new = zeros(n_clusters, n_data);
        for k = 1:n_clusters
            U_new(k, :) = 1 ./ (dist(k, :) ./ dist(1, :)) .^ 2;
        end
        U_new = U_new ./ repmat(sum(U_new, 1), n_clusters, 1);
        
        % 检查收敛
        if norm(U_new - U, 'fro') < tol
            U = U_new;
            break;
        end
        U = U_new;
    end
end