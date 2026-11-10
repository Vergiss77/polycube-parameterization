# Tasks

## 1. 高度图与简化算法（spec: polycube-simplification；design D1–D4）

- [x] 1.1 新建 `src/pcp/structure/height_simplify.{h,cpp}`，定义 `HeightGraph`、`SimplifyParams{enabled, max_stretch}`（含 to_json / update_from_json，未知键报错）和 `SimplifyResult`，并接入 CMake；验证：全量构建零错误，参数 JSON 往返单元测试通过
- [x] 1.2 实现单对合并代价：由相关边求可行区间 [a, b]，对 Σd 的导数做 100 次二分求 h*，代价取 Σd 的增量；验证：单元测试中两条边的解析最优高度误差 ≤ 1e-9，可行区间为空时报告不可行
- [x] 1.3 实现按轴的传递闭包（位集）、候选对枚举和贪心循环（同代价按组号排序），输出合并记录和每个面片的新高度；验证：单元测试覆盖链状 A → B → C 不合并 A、C、孤立组合并、不同轴不合并、全部边长度比位于 [1/s, s]、s 增大合并数不减、重复运行逐位一致
- [x] 1.4 实现 `build_height_graph(PatchGraph, heights)`：按逐位相等的高度和轴归组，由棱线 `lower_patch` / `upper_patch` 生成边，并统计缺端点的棱线；验证：单元测试中 L 形立方体组合体的组数、边数与手工计数一致，CLEVIS 的 `03_patch_graph.json` 构造出的初始切割面数等于每轴高度组数之和

- [x] 1.5 实现间隔边（design D1b）：在简化前的 polycube 上从每个三角形重心沿 ±轴向用 AABB 打射线，对第一个命中的同轴反向面片加"低 → 高"边并去重，统计间隔边数；验证：单元测试中 U 形 polycube 的底板底面和顶面、两个凸台内侧面之间都有间隔边且不被合并，只有两个顶面合并

- [x] 1.6 实现宽度边（design D1c）：在每个面片平面内，从三角形重心沿两条面内轴 ± 方向与面片边界求交，两侧边界面片同轴时加"低 → 高"边并去重；验证：单元测试中 L 形台阶面片两侧的两个同向面片之间有宽度边且不被合并

## 2. 按新高度重新映射（spec: Volumetric map with simplified heights；design D5）

- [x] 2.1 用 `heights_fixed` 模式串联"简化 → 重新映射"，计算简化后的切割面数、实测长度比范围、平面度误差和翻转数；验证：单元测试中 L 形组合体合并后面片节点的轴向坐标与新高度之差 ≤ 1e-9×对角线、翻转为 0、面片数与角点数不变、拓扑检查通过

- [x] 2.2 实现翻转回退（design D1d）：基准映射、全量检查、对合并前缀二分，统计贪心与保留的合并数和两种翻转数；验证：单元测试中只用棱线边建图的 U 形 polycube 只保留顶面合并、翻转为 0
- [x] 2.3 给最终映射的解缠加计算量上限 `untangle_budget`（design D1e，写入 `params.json`）；验证：单元测试中预算为 0 时解缠立即停止并如实报告剩余翻转数，默认预算下原有最终映射测试仍通过

## 3. 流水线与命令行接入（spec: Simplification stage in the pipeline / Simplification outputs；design D6）

- [x] 3.1 新增 `Stage::Simplify`（`--until simplify`，默认运行到此阶段）和 `PipelineParams.simplify`，写出 `05_simplify.json`、`05_polycube_simplified.tmsh/.vtu` 以及 `stats.json` 中的 `simplify` 字段；关闭时写 `{"enabled": false}`、不产生 `05_*`；验证：I-SPAR 上分别以默认、`--until structure`、关闭简化运行，检查输出文件有无、`params.json` 中的简化参数，以及顶层统计字段与改动前一致
- [x] 3.2 `pcp_cli` 新增 `--no-simplify`、`--max-stretch`，并更新 README（阶段列表、输出文件、新选项）；验证：`pcp_cli --help` 列出新选项，按 README 的命令在 I-SPAR 上运行，结果与描述一致

## 4. 模型验收与实验记录

- [x] 4.1 在 `model_ispar` 和 `model_clevis` 中增加简化检查：面片、角点、棱线数不变；切割面数不增（CLEVIS 严格减少）；长度比位于 [1/s, s]；平面度 ≤ 1e-9；翻转为 0；字段无 NaN；验证：`ctest -R '^model_(ispar|clevis)$'` 通过
- [x] 4.2 在 `model_tailsitter` 中增加简化检查：切割面数严格减少、翻转数 ≤ 0.01%；读取 `03_patch_graph.json`，分别以 s = 1.5 和 s = 3 直接调用简化，比较切割面数（较大 s 不多于较小 s）；验证：`ctest -R '^model_tailsitter$'` 通过
- [x] 4.3 在三个模型上对比 s = 1.5、2、3 的切割面数（每轴）、合并次数、简化后的平均与最大 ARAP、翻转数和耗时，确定默认 s，并在 `docs/experiments.md` 新增"高度简化"一节（含 D3 代价定义与论文的差异、TAILSITTER 台阶仍保留的说明、简化前后的对比图）；验证：文档中三个模型、三组 s 的数据完整
- [x] 4.4 在本地按 CI 命令序列复现（全新构建目录：构建、`ctest -L unit`、各模型测试）；推送后，CI 的构建加单元测试在两个平台通过，并手动触发一次带模型测试的运行；验证：本地全部通过，GitHub Actions 两个平台的任务通过

## Workflow follow-up

- 上一个 change `add-l1-polycube-pipeline` 还差 9.4（手动触发的 Windows 模型测试）。它通过后先归档上一个 change，再归档本 change，让 `polycube-simplification` 的主 spec 建在已有的主 spec 之上。
- 下一个 change：分块（沿简化后的高度组切割 polycube）、边界曲面拟合与三变量 B 样条体拟合（Chen §5.5.2 与 §6）。
