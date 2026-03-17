%% ------------------------------------------------------------------------------ %%
% apce.m by Jensentsts Wang
% aPCE模型函数库  v1.0.2
%
% 基本传入传出
% 以aPCE系数矩阵及其对应展开阶数表示一个aPCE模型
% coeffs: aPCE系数矩阵，行表示各阶展开，列表示基于alpha的各个基函数张量积结果
% x_matrix: 输入样本矩阵，行表示各样本，列表示各变量
% y_matrix: 输出样本矩阵，行表示各样本，列表示各输出变量
% degree: aPCE展开的阶数
% c_matrix: aPCE基函数系数矩阵，行表示各变量，列表示各阶展开的基函数系数。应当由pce_fit或pce_fit_auto的返回值保存以供后续预测使用。
% phi: aPCE基函数矩阵，行表示各样本，列表示各阶展开的基函数张量积结果
%
% 警告：
% ！！！建议传入归一化之后的值，避免运算过程中数字过大导致丢失精度！！！
% 归一化的方法：
% x_norm = apce().norm(x_matrix, x_min, x_max);
% 反归一化：
% y_fit = apce().abnorm(y_norm, y_min, y_max);
%
% 用法：
% [拟合阶数, apce拟合系数, apce基函数系数] = apce().fit_auto(最大允许的阶数, x_matrix, y_matrix, 允许的最大拟合误差)
% 如果无法在允许的最大拟合误差下完成拟合，返回的拟合系数为负数，也就是允许的最大拟合误差取负
% y_hat = apce().pred(拟合阶数, apce系数, apce基函数系数, x_fit)
% 如果训练数据为归一化的，那么输出也是归一化的。可以使用apce().norm()和apce().abnorm()进行相应的操作。

function F = apce
    %% 主要用到的函数：
    F.fit_auto = @apce_fit_auto;  % 自动拟合
    F.pred = @apce_predict;  % 预测
    F.norm = @apce_normalize;  % 归一化
    F.abnorm = @apce_abnormalize;  % 反归一化
    %% 提供更多的接口：
    F.fit = @apce_fit;  % 指定阶数的拟合
    F.tol = @apce_tol;  % 误差（L2）
    F.eval = @apce_eval;  % 给出此项，目的是用于效率优化
    %% 无特殊用途，用于检验目录下是否能找到此代码的文件
    F.test = @test;
end

function test()
    % 测试函数
    disp('This is a test function for aPCE module.');
end

%% ------------------------------------------------------------------------------ %%
% aPCE主函数 %

function val = apce_normalize(x_matrix, x_min, x_max)
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

function val = apce_abnormalize(x_matrix, x_min, x_max)
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

function [C,raw_mom] = apce_coef(po, x_vector)
    % input
    % po: polynomial order
    % x: data (1 column)
    %
    % output
    % C: coefficients of a polynomial (highest to lowest)
    % raw_mom: raw momoments used for getting C

    % k-th raw moments
    krm = zeros(1, 2 * po);
    for rm = 0 : 2 * po - 1
        krm(rm + 1) =  mean(x_vector.^rm);
    end

    % form matrices to get the coefficients
    A = zeros(po + 1, po + 1);
    for br = 1 : po
        A(br, :) = krm(br : br + po);
    end
    A(end, end) = 1;
    B = zeros(po + 1, 1);
    B(end) = 1;

    % determine coefficients
    C = flipud(A^(-1) * B);

    % output option
    if nargout == 2
        raw_mom = krm;
    end
end

function c_matrix = apce_gen_c(degree, x_matrix)
    n_samples = size(x_matrix, 1);
    n_vars = size(x_matrix, 2);
    c_matrix = zeros(n_vars, degree + 1, degree + 1);
    for deg_idx = 0 : degree
        for var_idx = 1 : n_vars
            c_matrix(var_idx, deg_idx + 1, 1:deg_idx + 1) = apce_coef(deg_idx, x_matrix(:, var_idx));
        end
    end
end

function phi_matrix = apce_gen_phi(degree, x_matrix, c_matrix)
    % 生成aPCE多变量的基函数矩阵
    % gen_basis_function: function handle，基函数, 参数为(degree, x_matrix, x_min, x_max)
    n_samples = size(x_matrix, 1);
    n_vars = size(x_matrix, 2);
    %% 预生成所有alpha
    n_alpha = 0;
    for d = 0 : degree
        n_alpha = n_alpha + gen_alpha_len(degree, n_vars);
    end
    alpha_list = zeros(n_alpha, n_vars);
    alpha_list_idx = 1;
    for d = 0 : degree
        alpha_d = gen_alpha(d, n_vars);
        alpha_list(alpha_list_idx : alpha_list_idx + size(alpha_d, 1) - 1, :) = alpha_d;
        alpha_list_idx = alpha_list_idx + size(alpha_d, 2);
    end
    %% 预生成所有基函数值
    phi_value = zeros(n_samples, n_vars, degree + 1);
    %{
    for deg_idx = 0 : degree
        phi_value(:, :, deg_idx + 1) = reshape(basis_function(deg_idx, x_matrix, x_min, x_max), size(x_matrix, 1), size(x_matrix, 2));
    end
    %}
    for deg_idx = 0 : degree
        for var_idx = 1 : n_vars
            c = reshape(c_matrix(var_idx, deg_idx + 1, 1:deg_idx + 1), 1, deg_idx + 1);
            phi = polyval(c, x_matrix(:, var_idx));
            phi_value(:, var_idx, deg_idx + 1) = phi;
        end
    end
    %% 求张量积
    phi_matrix = zeros(n_samples, n_alpha);
    % 注意phi_matrix的行表示各行表示样本，列表示各阶展开。
    for alpha_idx = 1 : n_alpha
        phi_matrix(:, alpha_idx) = basis_function_tensor_prod_from_alpha(alpha_idx);
    end

    %% 一些会用到的子函数
    function n_alpha = gen_alpha_len(d, n)
        % 计算n个变量、阶数之和<=p的所有多重指数alpha的数量
        n_alpha = nchoosek(d + n, n);
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
    function tp = basis_function_tensor_prod_from_alpha(alpha_idx)
        % 根据alpha生成基函数的张量积
        % 返回结果为一列向量，每行表示【对于该采样】的各变量的张量积，基于【此alpha列表】进行计算
        tp = ones(size(x_matrix, 1), 1);
        for idx_of_var = 1:size(x_matrix, 2)
            deg = alpha_list(alpha_idx, idx_of_var);
            tp = tp .* phi_value(:, idx_of_var, deg + 1);
        end
    end
end

function [coeffs, phi, c_matrix] = apce_fit(degree, x_matrix, y_matrix)
    % 计算aPCE拟合的系数
    % 将会返回phi矩阵，方便后续预测使用，提高效率
    c_matrix = apce_gen_c(degree, x_matrix);
    phi = apce_gen_phi(degree, x_matrix, c_matrix);
    coeffs = lsqminnorm(phi, y_matrix);
    %coeffs = double(lscov(phi, y_matrix));
end

function y_hat = apce_eval(degree, coeffs, c_matrix, phi, x_matrix)
    % 使用aPCE系数进行预测
    % 考虑到apce_fit_auto的效率优化，此处的phi、c_matrix可以为空

    if size(phi, 1) == 0
        c_matrix = apce_gen_c(degree, x_matrix);
        phi = apce_gen_phi(degree, x_matrix, c_matrix);
    end
    y_hat = phi * coeffs;
end

function y_hat = apce_predict(degree, coeffs, c_matrix, x_matrix)
    % 使用aPCE系数进行预测（不传入phi矩阵）
    % 一般是作为库以外的接口进行使用

    phi = apce_gen_phi(degree, x_matrix, c_matrix);
    y_hat = phi * coeffs;
end

function tol = apce_tol(degree, coeffs, c_matrix, phi, x_matrix, y_matrix)
    % 计算aPCE拟合的相对误差
    if size(phi, 1) ~= 0
        y_hat = phi * coeffs;
    else
        y_hat = apce_eval(degree, coeffs, c_matrix, [], x_matrix);
    end
    tol = norm(abs(y_matrix - y_hat)) / norm(y_matrix);
    %tol = max(abs(y_matrix - y_hat));
end

function [fitted_deg, coeffs, c_matrix] = apce_fit_auto(degree_max, x_matrix, y_matrix, tol_allow)
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
        [coeffs, phi, c_matrix] = apce_fit(degree, x_matrix, y_matrix);
        tol = apce_tol(degree, coeffs, c_matrix, phi, x_matrix, y_matrix);
        tol = max(tol);  % TODO
        fprintf("Tol at degree %d: %f\n", degree, tol);
    end
    %{
    if tol <= tol_allow
        fitted_deg = degree;
    else
        fitted_deg = -degree;
    end
    %}
    fitted_deg = degree;
end
