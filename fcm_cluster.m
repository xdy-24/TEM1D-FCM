 function [U, centers] = fcm_cluster(data, C, q, max_iter)
    % 输入：data [Nx2], C=聚类数, q=模糊指数
    % 输出：U [CxN]隶属度矩阵, centers [Cx2]聚类中心
    
    [N, dims] = size(data);
    if dims ~= 2
        error('输入数据的维度必须为Nx2');
    end
    
    U = rand(C, N);         % 初始化[CxN]
    U = U ./ sum(U, 1);     % 按列归一化
    
    for iter = 1:max_iter
        % 更新聚类中心
        centers = (U.^q * data) ./ sum(U.^q, 2);
        
        % 计算距离矩阵
        dist = pdist2(data, centers).^2;
        dist(dist == 0) = 1e-5;
        
        % 更新隶属度
        U_new = 1 ./ sum((dist ./ permute(dist, [3 2 1])).^(1/(q-1)), 3);
        U_new = squeeze(U_new)';
        U_new = U_new ./ sum(U_new, 1);  % 再次归一化
        
        % 收敛检查
        if max(abs(U_new(:) - U(:))) < 1e-5
            break;
        end
        U = U_new;
    end
end