function [is_polarized, cluster_type] = classify_layers(hpumct, centers, U, target_centers)
    % 分类层
    n_layers = size(hpumct, 1);
    is_polarized = false(n_layers, 1);
    cluster_type = zeros(n_layers, 1);
    
    % 参数归一化
    log_rho = log(hpumct(:,2));
    m_values = hpumct(:,4);
    
    % 判断每层的类型
    for i = 1:n_layers
        [max_u, cluster_idx] = max(U(:,i));
        cluster_type(i) = cluster_idx;
        
        % 判断是否为极化层
        is_polarized(i) = (cluster_idx == 2) && (max_u > 0.6);
    end
    
    % 输出分类结果到文件
    fid = fopen('layer_classification.txt', 'w');
    fprintf(fid, '层号\t厚度(m)\t电阻率(Ω·m)\t充电率\t隶属度\t极化层\n');
    fprintf(fid, '--------------------------------------------------------------\n');
    
    for i = 1:n_layers
        [max_u, cluster_idx] = max(U(:,i));
        fprintf(fid, '%d\t%.2f\t%.2f\t%.2f\t%.4f\t%s\n', ...
            i, hpumct(i,1), hpumct(i,2), hpumct(i,4), max_u, ...
            string(is_polarized(i)));
    end
    fclose(fid);
end 