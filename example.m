wtf;
%% PCE vs aPCE 拟合对比测试脚本
clear; clc; close all;

% 1. 构造测试数据 (模拟一个带有噪声的非线性系统)
% 假设函数: y = sin(pi * x1) + 0.5 * x2^2 + 随机噪声
n_samples = 100;
x_raw = rand(n_samples, 2) * 2 - 1; % 输入范围 [-1, 1]
y_raw = sin(pi * x_raw(:, 1)) + 0.5 * x_raw(:, 2).^2 + 0.05 * randn(n_samples, 1);

% 2. 准备测试集 (用于验证拟合效果)
n_test = 50;
x_test_raw = rand(n_test, 2) * 2 - 1;
y_test_real = sin(pi * x_test_raw(:, 1)) + 0.5 * x_test_raw(:, 2).^2;

% 数据归一化 (使用您库中的接口)
x_min = [-1, -1]; x_max = [1, 1];
x_norm = pce().norm(x_raw, x_min, x_max);
x_test_norm = pce().norm(x_test_raw, x_min, x_max);

%% 3. PCE 拟合
fprintf('--- 开始 PCE 自动拟合 ---\n');
max_deg = 5;
tol_req = 1e-2;
[pce_deg, pce_coeffs] = pce().fit_auto(max_deg, x_norm, y_raw, tol_req);
y_pred_pce = pce().pred(abs(pce_deg), pce_coeffs, x_test_norm);

%% 4. aPCE 拟合
fprintf('\n--- 开始 aPCE 自动拟合 ---\n');
[apce_deg, apce_coeffs, apce_c_mat] = apce().fit_auto(max_deg, x_norm, y_raw, tol_req);
y_pred_apce = apce().pred(abs(apce_deg), apce_coeffs, apce_c_mat, x_test_norm);

%% 5. 结果对比与可视化
% 计算误差 (RMSE)
pce_rmse = sqrt(mean((y_test_real - y_pred_pce).^2));
apce_rmse = sqrt(mean((y_test_real - y_pred_apce).^2));

fprintf('\n--- 结果总结 ---\n');
fprintf('PCE  最佳阶数: %d, 测试集 RMSE: %.4f\n', abs(pce_deg), pce_rmse);
fprintf('aPCE 最佳阶数: %d, 测试集 RMSE: %.4f\n', abs(apce_deg), apce_rmse);

% 绘图对比
figure('Color', 'w', 'Position', [100, 100, 1000, 400]);

% 子图1：PCE 预测 vs 真值
subplot(1, 2, 1);
plot(y_test_real, 'k-o', 'LineWidth', 1.5); hold on;
plot(y_pred_pce, 'r--x', 'LineWidth', 1.2);
title(['PCE 拟合 (Degree: ', num2str(abs(pce_deg)), ')']);
legend('真实值', 'PCE预测');
xlabel('测试样本编号'); ylabel('输出值');
grid on;

% 子图2：aPCE 预测 vs 真值
subplot(1, 2, 2);
plot(y_test_real, 'k-o', 'LineWidth', 1.5); hold on;
plot(y_pred_apce, 'b--d', 'LineWidth', 1.2);
title(['aPCE 拟合 (Degree: ', num2str(abs(apce_deg)), ')']);
legend('真实值', 'aPCE预测');
xlabel('测试样本编号'); ylabel('输出值');
grid on;

sgtitle('PCE 与 aPCE 代理模型性能对比');