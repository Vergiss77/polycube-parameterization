# Spec Delta

## Purpose

提供运行 polycube 流水线的命令行入口，统一管理参数、输出目录、中间结果与统计报告，便于在三个实验模型上复现与比较结果。

## ADDED Requirements

### Requirement: Run pipeline on a model
命令行工具 SHALL 接受一个输入表面网格路径与一个输出目录，依次执行校验、四面体化、ℓ1 变形与结构提取，并在成功时以零状态退出。

#### Scenario: Successful run
- **WHEN** 对 `models/CLEVIS-D64EEBFA.off` 指定输出目录运行工具
- **THEN** 工具以零状态退出，输出目录中包含各阶段的结果文件与统计报告

#### Scenario: Missing input
- **WHEN** 指定的输入文件不存在
- **THEN** 工具以非零状态退出并提示文件不存在

### Requirement: Stage selection and resume
工具 SHALL 支持只运行到指定阶段为止，以及从已保存的上一阶段结果继续运行，而无需重复计算前面的阶段。

#### Scenario: Stop after tetrahedralization
- **WHEN** 指定只运行到四面体化阶段
- **THEN** 输出目录中只有校验与四面体化结果，不执行变形

#### Scenario: Resume from saved tet mesh
- **WHEN** 指定从已保存的四面体网格开始运行变形
- **THEN** 工具跳过四面体化，结果与完整运行时一致

### Requirement: Configurable parameters
工具 SHALL 允许通过命令行选项或参数文件设置四面体化与变形的全部可调参数（如 β、α/ε 初值与倍率、最大阶段数、误差阈值、是否启用全局旋转），并把实际使用的完整参数写入输出目录。

#### Scenario: Parameters recorded
- **WHEN** 以自定义 β 运行工具
- **THEN** 输出目录中的参数记录包含该 β 值及其余所有参数的实际取值

### Requirement: Outputs for inspection
每次运行 SHALL 输出：带逐四面体 ARAP 畸变字段的 polycube 四面体网格（`.vtu`）、带面片标号的原始表面与 polycube 表面、面片图文件，以及机器可读的统计报告。

#### Scenario: Inspectable outputs
- **WHEN** 完整运行结束
- **THEN** 输出目录中存在上述五类文件，且 `.vtu` 文件可在 ParaView 中打开

### Requirement: Statistics report
统计报告 SHALL 至少包含：输入顶点/面数与亏格、四面体数与节点数、全局旋转、最终 polycube 误差、平均与最大 ARAP 畸变、面片数、角点数、翻转四面体数、阶段数、迭代次数与各阶段耗时。

#### Scenario: Report fields present
- **WHEN** 对 TAILSITTER 完整运行
- **THEN** 统计报告包含上述全部字段且数值合法（非 NaN）

### Requirement: Progress logging
工具 SHALL 在运行过程中输出每个调度阶段的 α、ε、polycube 误差与 ARAP 畸变，并支持安静模式与详细模式。

#### Scenario: Per-stage log
- **WHEN** 以默认日志级别运行变形
- **THEN** 每个调度阶段结束时输出一行包含 α、ε、polycube 误差与畸变的日志
