# Spec Delta

## Purpose

在浏览器中加载并交互查看流水线输出的三变量 B 样条多片体（G+Smo MultiPatch XML），显示离散后的边界面网格、缩放 Jacobian、内部接口面和控制网，并支持两个结果并排对比，用于检查和比较参数化结果。

## ADDED Requirements

### Requirement: Load a multipatch file
查看器 SHALL 能够加载流水线输出的 G+Smo MultiPatch XML（`07_bspline.xml`），加载方式包括：
- 文件选择对话框；
- 把文件拖进页面；
- 页面地址参数给出的、相对于静态服务根目录的路径。

文件内容不完整或不符合格式（缺少 MultiPatch、控制点数与节点向量不匹配、接口引用不存在的块等）时，查看器 MUST 在页面上显示说明原因的错误信息，不能崩溃，也不能显示部分结果。

#### Scenario: Load by file picker
- **WHEN** 用文件选择对话框打开 CLEVIS 的 `07_bspline.xml`
- **THEN** 页面显示 41 个块的边界面网格，信息面板显示块数 41、接口数 62、边界面数 122

#### Scenario: Load from URL
- **WHEN** 以仓库根目录为静态服务根，在页面地址中给出某个输出目录下 `07_bspline.xml` 的相对路径
- **THEN** 页面打开后直接显示该文件，不需要手动选择

#### Scenario: Invalid file
- **WHEN** 加载一个控制点数少于节点向量要求的文件
- **THEN** 页面显示错误信息，指出出错的块和原因，不显示任何几何

### Requirement: Evaluation consistent with the pipeline
查看器求值样条体得到的位置和 Jacobian SHALL 与 C++ 流水线一致：在同一组参数点上，位置之差不超过包围盒对角线的 1e-12 倍。信息面板中的缩放 Jacobian 统计 MUST 与流水线的统计报告一致：评估点数、非正点数、含非正点的块数相同，最小值与平均值之差不超过 1e-9。统计在与流水线相同的评估网格上计算（每块每方向为单元数 × 4 + 1 个点）。

#### Scenario: Positions match the C++ samples
- **GIVEN** C++ 对一个小算例写出的 `07_bspline.xml` 和 `07_bspline_samples.vtu`
- **WHEN** 用查看器的求值代码在相同的评估网格上求值
- **THEN** 每个点的位置与采样网格中对应点之差不超过包围盒对角线的 1e-12 倍，缩放 Jacobian 之差不超过 1e-9

#### Scenario: Statistics match the pipeline
- **WHEN** 加载 TAILSITTER 的 `07_bspline.xml`
- **THEN** 信息面板显示的评估点数为 67 750、非正点数为 29、含非正点的块数为 15，与 `stats.json` 的 `fit` 字段一致

### Requirement: Boundary surface mesh
查看器 SHALL 把文件中列出的每个边界面离散成三角网格显示，每个单元在每个方向上的细分数可调（1 到 32，默认 8）。相邻块在公共棱上的离散点 MUST 位置一致，不出现裂缝。可以单独开关网格线；网格线有两种：单元边界线，以及全部细分线。

#### Scenario: Change subdivision
- **WHEN** 把细分数从 8 改为 2
- **THEN** 显示的三角形数减少为原来的约 1/16，几何形状不变

#### Scenario: Element lines
- **WHEN** 打开单元边界线
- **THEN** 每个边界面上显示样条的单元划分（节点线）

### Requirement: Coloring and Jacobian markers
查看器 SHALL 提供三种着色方式：按块（相邻块颜色不同）、按缩放 Jacobian（连续色标，并显示色标和数值范围）、纯色。查看器 MUST 能把评估网格上缩放 Jacobian 非正的点标记出来（包括位于体内部的点），并可以单独开关这些标记。

#### Scenario: Jacobian coloring
- **WHEN** 对 TAILSITTER 选择按缩放 Jacobian 着色
- **THEN** 边界面按其上各点的缩放 Jacobian 着色，色标显示当前范围

#### Scenario: Non-positive markers
- **WHEN** 对 TAILSITTER 打开非正点标记
- **THEN** 显示 29 个标记点，位置与流水线评估网格中非正点的位置一致

### Requirement: Interior interfaces and clipping
查看器 SHALL 能够显示块与块之间的内部接口面（与边界面使用同样的离散方式，可单独开关、半透明显示），并提供沿 x、y、z 任一轴的剖切平面（位置可拖动调节，可翻转保留的一侧）。

#### Scenario: Show interfaces
- **WHEN** 对 CLEVIS 打开内部接口面
- **THEN** 62 个内部接口面被显示出来，可以看到块的划分

#### Scenario: Clip
- **WHEN** 沿 z 轴设置剖切平面并拖到模型中部
- **THEN** 平面一侧的几何被隐藏，配合内部接口面可以看到内部的块结构

### Requirement: Control net display
查看器 SHALL 能够显示控制点和控制网（同一块内沿三个参数方向相邻控制点之间的连线），可以显示全部块，也可以只显示选中的块。

#### Scenario: Control net of one block
- **WHEN** 点选一个块并打开"只显示选中块的控制网"
- **THEN** 只显示这个块的控制点和控制网线，控制点数为该块三个方向控制点数之积

### Requirement: Block picking
在三维视图中点击某个块的表面时，查看器 SHALL 高亮该块，并显示它的编号、三个方向的次数与单元数、控制点数，以及该块评估网格上缩放 Jacobian 的最小值和非正点数。

#### Scenario: Pick a block
- **WHEN** 点击 I-SPAR 的某个块
- **THEN** 该块被高亮，信息面板显示其编号和缩放 Jacobian 最小值（约为 1）

### Requirement: Side-by-side comparison
查看器 SHALL 支持同时加载两个文件，在左右两个视图中并排显示，两边各有自己的信息面板。默认两个视图的相机同步（旋转、平移、缩放一致），也可以取消同步。显示设置（细分数、着色方式、剖切位置等）MUST 同时作用于两个视图。

#### Scenario: Compare two results
- **WHEN** 左边加载默认参数的 CLEVIS 结果，右边加载以 `--element-size 0.025` 得到的结果
- **THEN** 两个视图并排显示，旋转一边另一边同步转动，两个信息面板分别显示各自的控制点数

### Requirement: Offline operation
查看器 SHALL 在没有网络连接时正常工作：所有脚本和样式都来自仓库内的文件，运行时不访问任何外部地址。引入的第三方代码 MUST 使用固定版本，并附带其许可证文件。

#### Scenario: No network
- **WHEN** 断开网络后，用本地静态服务打开查看器并加载文件
- **THEN** 页面功能完全正常，浏览器开发者工具中没有对外部地址的请求

### Requirement: Usable on the provided models
在开发机上以默认细分数加载 TAILSITTER 的结果时，查看器 SHALL 在 10 秒内完成显示，之后旋转等交互保持流畅（不出现明显卡顿）。

#### Scenario: Load TAILSITTER
- **WHEN** 在开发机上加载 TAILSITTER 的 `07_bspline.xml`
- **THEN** 10 秒内显示完成，信息面板给出加载与离散的耗时
