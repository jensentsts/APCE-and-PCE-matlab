# PCE 与 aPCE 工具箱 (MATLAB)
[English](./)
本项目提供了**多项式混沌展开 (Polynomial Chaos Expansion, PCE)** 与**任意多项式混沌展开 (arbitrary Polynomial Chaos Expansion, aPCE)** 的 MATLAB 实现。  
通过给定的输入‑输出样本，可以自动或手动拟合 PCE / aPCE 模型，并用于预测、误差分析等任务。

---

## 目录
- [简介](#简介)
- [文件说明](#文件说明)
- [安装](#安装)
- [快速开始](#快速开始)
  - [1. 数据归一化](#1-数据归一化)
  - [2. 拟合 PCE 模型](#2-拟合-pce-模型)
  - [3. 预测](#3-预测)
  - [4. 自动选择阶数](#4-自动选择阶数)
- [函数文档](#函数文档)
  - [pce.m 函数](#pcem-函数)
  - [apce.m 函数](#apcem-函数)
- [注意事项](#注意事项)
- [贡献](#贡献)
- [许可证](#许可证)
- [作者](#作者)

---

## 简介

**多项式混沌展开 (PCE)** 是一种将模型输出表示为输入随机变量正交多项式级数展开的方法，广泛应用于不确定性量化、灵敏度分析等领域。  
**任意多项式混沌展开 (aPCE)** 是 PCE 的推广，它利用输入数据的统计矩构造与数据分布自适应匹配的正交多项式基，无需假设输入服从特定分布（如均匀、正态）。

本工具箱提供：
- 基于 Legendre 多项式的标准 PCE（输入需归一化至 `[-1, 1]`）
- 基于数据矩构造正交多项式的 aPCE（可处理任意分布输入）
- 自动选择展开阶数的功能
- 归一化 / 反归一化辅助函数
- 误差（L2 相对误差）计算

---

## 文件说明

| 文件名    | 描述                                           |
| --------- | ---------------------------------------------- |
| `pce.m`   | 标准 PCE 实现（基函数：Legendre 多项式）       |
| `apce.m`  | 任意 PCE 实现（基函数由输入数据矩构造）        |

两个文件的接口风格一致，方便切换使用。

---

## 安装

1. 将 `pce.m` 和 `apce.m` 下载到本地任意文件夹。
2. 在 MATLAB 中，将该文件夹添加到搜索路径：
   ```matlab
   addpath('文件夹路径');
   ```
   或者右键文件夹 → **添加到路径** → **选定的文件夹**。

---

## 快速开始

以下示例演示如何使用 aPCE 拟合一个简单的一维函数 `y = x^2`，并预测新点。

### 1. 数据归一化
输入通常需要归一化到 `[-1, 1]` 以避免数值问题。`apce()` 提供了 `norm` 和 `abnorm` 函数。

```matlab
% 生成训练数据
x = rand(100, 1) * 4 - 2;        % x ∈ [-2, 2]
y = x.^2;

% 归一化输入
x_min = min(x);
x_max = max(x);
x_norm = apce().norm(x, x_min, x_max);
```

### 2. 拟合 aPCE 模型
使用 `fit_auto` 自动选择合适阶数（最大允许阶数设为 10，允许误差 1e-3）：

```matlab
[degree, coeffs, c_matrix] = apce().fit_auto(10, x_norm, y, 1e-3);
```
- `degree`：最终采用的展开阶数
- `coeffs`：PCE 系数矩阵（每一列对应一个输出）
- `c_matrix`：aPCE 特有的基函数系数矩阵（用于预测）

### 3. 预测
对新点进行预测（需先归一化新点）：

```matlab
x_new = [-1.5; 0; 1.8];
x_new_norm = apce().norm(x_new, x_min, x_max);
y_pred_norm = apce().pred(degree, coeffs, c_matrix, x_new_norm);
y_pred = apce().abnorm(y_pred_norm, min(y), max(y));   % 反归一化输出
```

### 4. 自动选择阶数
`fit_auto` 会从阶数 0 开始递增，直到拟合误差低于 `tol_allow` 或达到 `degree_max`。若成功，返回正阶数；若失败，返回负值（绝对值表示尝试的最大阶数）。

---

## 函数文档

### pce.m 函数

所有函数通过 `pce()` 返回的结构体调用，例如 `pce().fit_auto(...)`。

| 函数名 | 说明 |
|--------|------|
| `fit_auto(degree_max, x, y, tol_allow)` | 自动选择阶数并拟合，返回 `[fitted_deg, coeffs]` |
| `pred(degree, coeffs, x_new)` | 用训练好的 PCE 模型预测新点 |
| `norm(x, x_min, x_max)` | 将数据线性映射到 `[-1, 1]`，若省略 `x_min, x_max` 则使用数据本身的最值 |
| `abnorm(x, x_min, x_max)` | 将 `[-1, 1]` 数据映射回原始区间 |
| `legendre(n, x)` | 计算 n 阶 Legendre 多项式在 x 处的值（自动处理归一化） |
| `fit(degree, x, y)` | 用指定阶数拟合 PCE，返回 `[coeffs, phi]` |
| `tol(degree, coeffs, phi, x, y)` | 计算相对 L2 误差 |
| `eval(degree, coeffs, phi, x)` | 用已有 phi 矩阵预测（效率优化） |
| `test()` | 测试函数，仅显示信息 |

### apce.m 函数

通过 `apce()` 结构体调用。

| 函数名 | 说明 |
|--------|------|
| `fit_auto(degree_max, x, y, tol_allow)` | 自动选择阶数并拟合，返回 `[fitted_deg, coeffs, c_matrix]` |
| `pred(degree, coeffs, c_matrix, x_new)` | 用训练好的 aPCE 模型预测新点 |
| `norm(x, x_min, x_max)` | 同 pce |
| `abnorm(x, x_min, x_max)` | 同 pce |
| `fit(degree, x, y)` | 用指定阶数拟合 aPCE，返回 `[coeffs, phi, c_matrix]` |
| `tol(degree, coeffs, c_matrix, phi, x, y)` | 计算相对 L2 误差 |
| `eval(degree, coeffs, c_matrix, phi, x)` | 用已有 phi 矩阵预测（效率优化） |
| `test()` | 测试函数 |

---

## 注意事项

1. **强烈建议对输入数据进行归一化**（映射到 `[-1, 1]`），否则高次项可能导致数值溢出或严重精度损失。工具箱提供了 `norm` / `abnorm` 辅助函数。
2. aPCE 的基函数由输入数据的一至二阶矩构造，因此训练数据应能代表输入的真实分布。
3. 当输入维度较高或阶数较大时，基函数数量会急剧增加（组合爆炸），请合理设置最大阶数。
4. `fit_auto` 中的 `tol_allow` 为相对 L2 误差阈值，可根据实际问题调整。
5. 代码使用 `lsqminnorm` 求解最小二乘问题，若出现警告可考虑使用 `lscov` 替代（注释中已提供备选）。

---

## 贡献

欢迎通过 Issues 或 Pull Requests 提交改进建议、bug 修复或新功能。请确保代码风格与现有文件一致，并添加必要的注释。

---

## 许可证

本项目采用 MIT 许可证。详情请参阅 [LICENSE](LICENSE) 文件。

---

## 作者

Jensentsts Wang  
如有问题，可联系：jensentsts.wang@example.com（请替换为真实邮箱）

---

*最后更新：2026-03-18*
