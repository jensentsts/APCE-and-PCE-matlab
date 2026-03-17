%% ------------------------------------------------------------------------------ %%
% pce.m by Jensentsts Wang
% PCE模型函数库  v1.3.0
%
% 基本传入传出
% 以PCE系数矩阵coeffs及其对应展开阶数degree表示一个PCE模型
% coeffs: PCE系数矩阵，行表示各阶展开，列表示基于alpha的各个基函数张量积结果
% x_matrix: 输入样本矩阵，行表示各样本，列表示各变量
% y_matrix: 输出样本矩阵，行表示各样本，列表示各输出变量
% degree: PCE展开的阶数
% phi: PCE基函数矩阵，行表示各样本，列表示各阶展开的基函数张量积结果
%
% ！！！注意：必须传入归一化之后的值，避免运算过程中数字过大导致丢失精度
% 归一化的方法：
% x_norm = pce().norm(x_matrix, x_min, x_max);
% 反归一化：
% y_fit = pce().abnorm(y_norm, y_min, y_max);
%
% [拟合阶数, pce拟合系数] = pce().fit_auto(最大允许的阶数, x_matrix, y_matrix, 允许的最大拟合误差)
% 如果无法在允许的最大拟合误差下完成拟合，返回的拟合系数为负数，也就是允许的最大拟合误差取负
% y_hat = pce().pred(拟合阶数, pce系数, x_fit)
% 如果训练数据为归一化的，那么输出也是归一化的。可以使用pce().norm()和pce().abnorm()进行相应的操作。

function F = pce
    % 主要用这两个：
    F.fit_auto = @pce_fit_auto;
    F.pred = @pce_predict;

    F.norm = @pce_normalize;
    F.abnorm = @pce_abnormalize;
    F.legendre = @pce_legendre;
    F.fit = @pce_fit;
    F.tol = @pce_tol;
    F.eval = @pce_eval;

    F.test = @test;
end

function test()
    % 测试函数
    disp('This is a test function for PCE module.');
end

%------------------------------------------------------------------------------%
% PCE主函数 %
%------------------------------------------------------------------------------%

function val = pce_normalize(x_matrix, x_min, x_max)
    % 归一化函数(线性映射到[-1, 1])
    
    if length(x_min) ~= length(x_max)
        error("x_min与x_max长度不一致！");
    end

    % 支持单输入
    if nargin == 1
        x_min = zeros(size(x_matrix, 2));
        x_max = ones(size(x_matrix, 2)) * inf;
        for i = 1:size(x_matrix, 2)
            x_min(i) = min(x_matrix(:, i));
            x_max(i) = max(x_matrix(:, i));
        end
    end

    val = x_matrix;
    if length(x_min) > 1
        for i = 1:length(x_min)
            val(:, i) = (x_matrix(:, i) - x_min(i)) / (x_max(i) - x_min(i)) * 2 - 1;
        end
    else
        val = (x_matrix - x_min) / (x_max - x_min) * 2 - 1;
    end
end

function val = pce_abnormalize(x_matrix, x_min, x_max)
    % 反归一化函数(线性映射回原区间)
    
    if length(x_min) ~= length(x_max)
        error("x_min与x_max长度不一致！");
    end

    % 支持单输入
    if nargin == 1
        x_min = zeros(size(x_matrix, 2));
        x_max = ones(size(x_matrix, 2)) * inf;
        for i = 1:size(x_matrix, 2)
            x_min(i) = min(x_matrix(:, i));
            x_max(i) = max(x_matrix(:, i));
        end
    end

    val = x_matrix;
    if length(x_min) > 1
        for i = 1:length(x_min)
            val(:, i) = ( (x_matrix(:, i) + 1) / 2 ) * (x_max(i) - x_min(i)) + x_min(i);
        end
    else
        val = ( (x_matrix + 1) / 2 ) * (x_max - x_min) + x_min;
    end
end


function val = pce_legendre(n, x_matrix)
    % 计算Legendre多项式、自动归一化
    leg_vals = legendre(n, x_matrix);
    if n == 0
        val = leg_vals;
    else
        val = leg_vals(1, :);
    end
end

function phi_matrix = pce_gen_phi(degree, bf, x_matrix)
    % 生成PCE多变量的基函数矩阵
    % bf: function handle，基函数, 参数为(degree, x_matrix)

    % 效率优化：预生成所有alpha
    alpha_len = 0;
    for d = 0 : degree
        alpha_len = alpha_len + gen_alpha_len(degree, size(x_matrix, 2));
    end
    alpha_list = zeros(alpha_len, size(x_matrix, 2));
    alpha_list_idx = 1;
    for d = 0 : degree
        alpha_d = gen_alpha(d, size(x_matrix, 2));
        alpha_list(alpha_list_idx : alpha_list_idx + size(alpha_d, 1) - 1, :) = alpha_d;
        alpha_list_idx = alpha_list_idx + size(alpha_d, 2);
    end

    % 效率优化：预生成所有基函数值
    phi_value = zeros(size(x_matrix, 1), size(x_matrix, 2), degree + 1);
    for deg_idx = 0 : degree
        phi_value(:, :, deg_idx + 1) = reshape(bf(deg_idx, x_matrix), size(x_matrix, 1), size(x_matrix, 2));
    end

    phi_matrix = zeros(size(x_matrix, 1), alpha_len);
    % 注意phi_matrix的行表示各行表示样本，列表示各阶展开。
    for alpha_idx = 1 : alpha_len
        phi_matrix(:, alpha_idx) = gen_bf_tensor_prod_from_alpha(alpha_idx);
    end

    function alpha_len = gen_alpha_len(d, n)
        % 计算n个变量、阶数之和<=p的所有多重指数alpha的数量
        alpha_len = nchoosek(d + n, n);
    end

    function alpha = gen_alpha(d, n)
        % 生成n个变量、阶数之和<=p的所有多重指数alpha
        % 例如： p=3, n=2 时，生成的alpha为：
        % 0 0
        % 0 1
        % 0 2
        % 0 3
        % 1 0
        % 1 1
        % ...
        % 3 0
        % 例如2 2，违反规则，不会被生成
        alpha = [];
        current = zeros(1, n);
        recursive_gen(current, 1, d);

        function recursive_gen(current, i, remaining)
            % 一种递归生成alpha的算法
            if i > n
                alpha = [alpha; current];
                return;
            end
            for ai = 0:remaining
                current(i) = ai;
                recursive_gen(current, i + 1, remaining - ai);
            end
        end
    end

    function tp = gen_bf_tensor_prod_from_alpha(alpha_idx)
        % 根据alpha生成基函数的张量积
        % 返回结果为一列向量，每行表示【对于该采样】的各变量的张量积，基于【此alpha列表】进行计算
        tp = ones(size(x_matrix, 1), 1);
        for var_idx = 1:size(x_matrix, 2)
            deg = alpha_list(alpha_idx, var_idx);
            tp = tp .* phi_value(:, var_idx, deg + 1);
        end
    end
end

function [coeffs, phi] = pce_fit(degree, x_matrix, y_matrix)
    % 计算PCE拟合的系数
    % 将会返回phi矩阵，方便后续预测使用，提高效率

    % 此处直接使用了pce_legendre作为基函数，后续也可以替换为其他基函数，算作一个TODO了
    phi = pce_gen_phi(degree, @pce_legendre, x_matrix);
    coeffs = lsqminnorm(phi, y_matrix);
end

function y_hat = pce_eval(degree, coeffs, phi, x_matrix)
    % 使用PCE系数进行预测
    if size(phi, 1) == 0
        phi = pce_gen_phi(degree, @pce_legendre, x_matrix);
    end
    y_hat = phi * coeffs;
end

function y_hat = pce_predict(degree, coeffs, x_matrix)
    % 使用PCE系数进行预测（不传入phi矩阵）
    phi = pce_gen_phi(degree, @pce_legendre, x_matrix);
    y_hat = phi * coeffs;
end

function tol = pce_tol(degree, coeffs, phi, x_matrix, y_matrix)
    % 计算PCE拟合的相对误差
    if size(phi, 1) ~= 0
        y_hat = phi * coeffs;
    else
        y_hat = pce_eval(degree, coeffs, [], x_matrix);
    end
    tol = norm(abs(y_matrix - y_hat)) / norm(y_matrix);
    %tol = max(abs(y_matrix - y_hat));
    fprintf("Tol at degree %d: %f\n", degree, tol);
end

function [fitted_deg, coeffs] = pce_fit_auto(degree_max, x_matrix, y_matrix, tol_allow)
    % 自动确定合适的拟合阶数
    % fitted_deg: 最终拟合的阶数。返回正数表示成功拟合，负数表示未能在最大阶数内达到要求。为此，fitted_deg的绝对值表示尝试的阶数或达到最大阶数。

    % 防呆设计，参数过滤
    if degree_max == inf
        degree_max = 1e3;
    end

    tol = inf;
    degree = -1;
    while tol > tol_allow && degree <= degree_max
        degree = degree + 1;
        [coeffs, phi] = pce_fit(degree, x_matrix, y_matrix);
        tol = pce_tol(degree, coeffs, phi, x_matrix, y_matrix);
    end
    if tol <= tol_allow
        fitted_deg = degree;
    else
        fitted_deg = -degree;
    end
end
