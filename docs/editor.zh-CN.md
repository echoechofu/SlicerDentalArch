# DentalArchEditor：Slicer 牙弓与拥挤度修订模块

[English](editor.md) · [项目首页](../README.zh-CN.md)

本模块用于在 Slicer 中人工修订近远中冠宽测量点和参考牙弓，旁边同步显示空间分析结果。
算法初始值与当前工作副本分开保存，修订结果另存为新版本。

这是测量编辑原型，输入为已经标注的场景和测量元数据；它不会从任意 STL 自动识别或分割牙齿。
尚未验证临床诊断准确度。所有长度使用输入模型的原始坐标单位，不自动认定为毫米。
源码采用 BSD-3-Clause，尚未发布到 Slicer 扩展商店。公开源码包不包含病例模型或病例坐标。

## 本地安装

目标环境为 **3D Slicer 5.12.4**，使用其自带的 NumPy、SciPy、VTK 及 Qt/PythonQt。
普通系统 Python 可以运行数值测试，但不能单独启动本模块的 Slicer 界面。

1. 将源码包解压到一个固定位置。
2. 打开 Slicer 的 **Edit → Application Settings → Modules**。
3. 在 **Additional module paths** 中加入包含 `DentalArchEditor.py` 的目录，
   即源码包中的 `DentalArchEditor` 子目录，不是源码包最外层目录。
4. 重启 Slicer，从模块选择器中选择 **Dentistry → Dental Arch Editor / 牙弓修订**。

这种 Python 模块安装方式无需编译 Slicer。参见
[Slicer 官方扩展说明](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html)。

## 先试合成示例

点击 **加载合成示例 / Synthetic demo**，模块生成数学曲线和简化球形牙冠。
可以用它练习拖点、调弓、撤销／重做、查看实时数值及保存修订版，无需提供真实病例。
这些形状和测量值是人为构造的，不代表解剖示例。

## 打开新流程生成的场景

[算法流程](algorithms.zh-CN.md)生成的 `DentalArchScene.mrb` 已经初始化编辑器状态。
打开该 MRB 并进入本模块即可复核，无需再次导入 JSON。MRB 保存数据，不保存运行中的
Python 回调，因此重新加载后，需要安装并打开本模块来恢复实时计算。

## 导入旧版标注场景

先打开已经标注的 Slicer 场景，再使用 **导入原算法数据 / Import baseline**，选择对应的两份 JSON：

- 近远中结果 JSON：牙齿标识、M/D 原始坐标、坐标框架及源模型标识。
- 参考牙弓结果 JSON：原始 B 样条控制点、节点向量、阶数、参考平面及测量范围。

第一版导入器兼容本项目的测量结果格式。其他软件需要输出相同格式或编写适配器。
只有一个 STL 或一串没有牙齿身份的点，无法确定完整测量协议。
现有场景中的 501 个弓线显示采样点不会变成 501 个编辑手柄；模块读取真实的样条参数，
初次导入时保留原曲线，不暗中重拟合。

通过牙冠选择器，把当前牙齿关联到其独立牙冠模型，再启用表面吸附。
`VL/VR` 沿用输入数据的命名约定，模块不会据此推断患者解剖左右或 FDI 编号。

其他流程编写适配器时，主要字段如下：

| 文件 | 测量字段 |
| --- | --- |
| 近远中 JSON | `source_sha256`、`pca_center`、`pca_axes_columns`；`teeth` 列表中含 `id`、`mesial_native`、`distal_native`，可附 `kind` |
| 参考弓 JSON | `source_STL_sha256`、相同坐标框架、`reference_plane_aligned_z`、`control_points_xy`、`knots`、`degree`、`input_tooth_ids`、`endpoint_ids`、`projected_arch_length` |

两份文件的源标识须一致。若给出 `input_v2_json_sha256`，它必须匹配近远中 JSON 的实际文件字节。
端点名称格式为 `牙齿ID_M` 或 `牙齿ID_D`。框架矩阵按列存储右手正交轴，
局部坐标等于 `(点 - pca_center) @ pca_axes_columns`。
导入器识别 `牙齿ID_MD`／`牙齿ID_MD_CANDIDATE` 基线线节点及
`牙齿ID_crown`／`牙齿ID_isolated_crown` 牙冠模型，同名歧义会拒绝导入。
修订版保存后使用持久节点引用恢复关系，不依赖用户界面中的显示名称。

## 人工修订

1. 选择牙齿，进入 **调整接触测量点 / Edit M-D points**，直接拖动三维视图中的 M 或 D 点。
   启用吸附时，点投影到该牙关联的牙冠。两点当前距离作为修订后的牙宽。
2. 进入 **调整牙弓手柄 / Edit arch handles**，在固定参考平面内拖动内部弓形控制柄。
   B 样条控制柄不一定在实际弓线上。
3. 侧栏显示所需空间 `R`、可用空间 `A`、有符号拥挤量 `C=R-A`，并对照原始值、当前值和变化量。
4. 使用 **撤销 / Undo、重做 / Redo、恢复当前牙原始点位 / Reset tooth** 和原始对比显示复核修改。
5. 需要依据当前接触测量点重新生成参考弓时，主动点击 **Refit arch／重新拟合牙弓**。

默认手动修订中，改变内部牙齿的 M/D 点会更新 `R`，内部弓形控制点保持不变。
两侧第一磨牙近中点仍定义弓线边界：调整它们时，相应投影端点移动，`A` 也会更新。
编辑期间参考坐标框架保持固定。重新拟合是单独操作，避免改一个测量点时替换已经人工调好的弓形。

当前支持平面三次 B 样条，并逐段检查其沿两边界连线方向的推进保持单调。
这能排除自交，也会保守地拒绝一些没有自交但横向回折的曲线。反向弯曲作为复核提示。
输入无效时清空当前合计值，不把上一次结果冒充为当前结果。
标注和牙冠应保留在共同的原始坐标系中，不通过父级变换来另行改变测量框架。

自动初始化的后牙上 1/3 和方向规则用于复核提示。手动点不会被暗中拉回算法的最大截面弦位置。
接受修订前，应查看点位规则偏差及几何状态提示。

## 数学定义

默认上颌测量范围为两侧第一磨牙近中点之间：

```text
R = 两侧第二前磨牙之间共 10 颗牙的 ||M_i-D_i|| 之和
A = 两个第一磨牙近中投影端点之间的平面参考样条弧长
C = R-A：正值表示净空间不足，负值表示净空间余量
D = A-R：正值表示净空间余量，负值表示净空间不足
```

第一磨牙只设边界，其宽度不计入 `R`；第二磨牙不属于此默认范围。
总量略有余隙不代表没有局部旋转、重叠或错位。近远中点是牙冠宽度的测量位置，
不保证是现状下邻牙之间的真实接触点。参考弓是测量基准，不是理想排牙或治疗目标。

## 保存和恢复

使用模块的 **Save new revision／保存为新修订版**，选择目录。
每次修订的 MRB 场景、JSON 状态和 CSV 测量表应保存在一起。
该操作拒绝覆盖已有修订文件或导入的原始文件。

修订版保留基线与当前点位、真实样条控制点、参考平面、计算范围及结果。
MRB 内包含查看所需的模型。安装模块后重新打开修订场景，即可继续交互计算。
MRB 保存的是数据，不会保存运行中的 Python 监听函数；模块在重新加载时恢复监听关系。

模块的防覆盖规则适用于它自己的修订保存按钮。Slicer 通用 **Save** 窗口是独立功能，
不要在其中主动指定原始文件作为覆盖目标。
CSV 是保存时的计算快照，在外部修改 CSV 不会同步修改场景。

## 开发和测试

纯数值代码独立于 Slicer，测试只使用合成数据：

```sh
python -m pip install numpy scipy vtk
python -m unittest discover -s tests -p "test_*.py" -v
```

GitHub Actions 执行这一数值测试集。通过数值测试不等于交互或临床准确度验证。
还需在 Slicer 中检查拖动、所属牙冠吸附、撤销／重做、端点联动、重新拟合、场景保存重载和基线保护。

源码另附不使用病例数据的场景集成测试：

```sh
Slicer --no-main-window --no-splash --python-script DentalArchEditor/Testing/Python/SceneRoundTrip.py
```

它检查编辑、撤销、另存新文件、原始值保留及重载后的自动监听恢复，仅生成临时合成数据。
当前版本的本地验证范围见 [VALIDATION.md](../VALIDATION.md)。

CMake 文件提供 Slicer 扩展构建入口；本地优先使用附加模块路径安装。
打包和跨平台兼容性需在目标环境分别验证。准备提交 Slicer Extensions Index 时，
再核对扩展元数据并补充所需的公开图标和文档地址，按
[官方发布流程](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html)办理。

贡献说明见 [CONTRIBUTING.md](../CONTRIBUTING.md)，依赖许可与方法背景见
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md)，项目许可见 [LICENSE](../LICENSE)。

