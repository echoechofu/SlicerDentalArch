# SlicerDentalArch

[English](README.md)

把**牙齿标注算法、Slicer 场景和交互修改模块**放在同一个项目中：算法生成牙尖、
磨牙中央窝、四颗切牙切缘、近远中测量点和参考牙弓；随后在 DentalArchEditor 中
调整测量点或参考牙弓，实时查看拥挤度，并另存修改结果，保留最初的算法基线。

```text
STL → 检查网格组件，人工确认牙位和坐标方向
    → 牙齿标注与参考牙弓算法 → DentalArchScene.mrb
    → DentalArchEditor 手动调整 → 新版本 MRB + JSON + CSV
```

当前是实验性的半自动流程：输入 STL 必须已有可分离的单牙牙冠网格组件，不能直接
对任意连成一体的口扫网格自动分牙。牙位和方向需要人工复核；软件测试不代表解剖
定位精度或临床准确性验证。所有长度保留 STL 原始坐标单位，STL 本身不提供毫米依据。

项目源码采用 BSD-3-Clause 许可，不包含患者数据。当前提供源码安装方式，尚未上架
Slicer Extensions Index。

## 运行条件

- 目标环境：**3D Slicer 5.12.4**，已在 Windows 本地验证。
- 使用 Slicer 自带的 Python、NumPy、SciPy、VTK 和 Qt/PythonQt 环境。
- 二进制 STL，包含 14 个可独立识别的上颌牙冠组件，可以同时含有其他底座组件。
  被选中的牙冠必须封闭、非退化，并保持一致的朝外面片方向。详见[算法说明](docs/algorithms.zh-CN.md)。
- 检测参数按牙齿在毫米量级坐标中的大小设置，包含固定长度的网格步长、平滑尺度和
  搜索阈值。遇到其他坐标尺度，需要先检查模型范围和参数；程序不会自动换算单位。

## 从 STL 到可修改的 Slicer 场景

### 1. 安装编辑模块

下载或克隆本仓库。在 Slicer 中打开 **Edit → Application Settings → Modules →
Additional module paths**，添加仓库内的 `DentalArchEditor` 文件夹，重启 Slicer，
然后选择 **Dentistry → Dental Arch Editor / 牙弓修订**。

[编辑器说明](docs/editor.zh-CN.md)介绍了面板操作、修订保存和已有标注场景的导入方式。

### 2. 检查模型，准备人工复核配置

使用 Slicer 安装目录中的 `PythonSlicer` 运行。以下为在仓库根目录执行的 PowerShell
示例，请替换软件路径和输入文件路径：

```powershell
& "C:\path\to\Slicer\bin\PythonSlicer.exe" algorithms/inspect_stl.py --input model.stl --output-dir work/inspect
```

检查组件统计和带组件编号的俯视、仰视图。把生成的
`work/inspect/config.template.json` 复制为私有的 `reviewed.json`，填写 14 个牙位对应的
组件编号，并确认坐标框架。只有完成检查后，才把两个复核标志设为 `true`。
PCA 只提供初始方向建议，不能替代对上下、前后和左右的解剖判断。

完整配置、适用条件和故障排查见[算法使用说明](docs/algorithms.zh-CN.md)。患者输入、
配置和结果请保存在公共源码目录之外。

### 3. 运行算法并生成场景

```powershell
& "C:\path\to\Slicer\Slicer.exe" --ignore-slicerrc --disable-settings --no-main-window --no-splash --python-script algorithms/run_in_slicer.py -- --input model.stl --config reviewed.json --output-dir output/session
```

主要产物是 **`output/session/DentalArchScene.mrb`**，其中封装牙冠模型、标注、参考牙弓、
算法基线和可编辑测量状态。同时输出 JSON/CSV 中间结果及运行清单。每次运行请使用
新的输出目录。可选 `--skip-sensitivity` 跳过切缘、近远中测量和参考牙弓的额外敏感性分析。

如果只需要 JSON/CSV，不创建场景：

```powershell
& "C:\path\to\Slicer\bin\PythonSlicer.exe" algorithms/run_all.py --input model.stl --config reviewed.json --output-dir output/session
```

### 4. 复核、调整和另存

打开生成的 MRB，进入 DentalArchEditor。场景已初始化，无需再次导入算法 JSON。
拖动牙齿近远中点或参考牙弓控制点，面板会更新所需长度、可用长度和两者差值。
调整内部牙齿测量点不会自动重拟合已修改的参考牙弓。

使用模块内的 **Save new revision**，在新修订目录中保存 MRB、JSON 和 CSV，保留基线与
当前结果。安装模块后重新打开修订场景，即可恢复实时计算。Slicer 自带的普通 Save
操作是另一个功能，仍可能覆盖用户指定的同名文件。

## 测量定义

默认范围为两侧第一磨牙近中点在参考平面上的投影之间，所需长度合计
**第二前磨牙到第二前磨牙，共 10 颗牙**。

| 量 | 定义 |
| --- | --- |
| 所需长度 `R` | 10 颗牙各自三维近远中测量线段长度之和 |
| 可用长度 `A` | 平面三次 B 样条参考牙弓的解析曲线弧长 |
| 拥挤度 `C = R − A` | 正值表示长度不足，负值表示净剩余间隙 |
| 间隙差 `D = A − R` | 与上述拥挤度相反的符号约定 |

参考牙弓按当前总体形态平滑局部错位，并非治疗目标牙弓。第一磨牙宽度不计入所需
长度，第二磨牙在测量范围之外。净剩余间隙不排除局部拥挤或旋转；程序不自动给出
临床严重程度分级。

## 示例与测试

模块中的 **加载合成示例 / Synthetic demo** 使用数学曲线和简单几何体生成演示场景，
无需患者数据。它用于验证编辑与保存流程，不用于证明牙齿标注精度。

```powershell
python -m pip install numpy scipy vtk
python -m unittest discover -s tests -v
& "C:\path\to\Slicer\Slicer.exe" --no-main-window --no-splash --python-script DentalArchEditor/Testing/Python/SceneRoundTrip.py
```

更多信息：[验证范围](VALIDATION.md)、[参与开发](CONTRIBUTING.md)、
[第三方依赖与方法参考](THIRD_PARTY_NOTICES.md)、[开源许可](LICENSE)。

| 目录 | 内容 |
| --- | --- |
| `algorithms/` | 输入检查、人工复核配置、标注/测量/拟合算法及场景生成 |
| `DentalArchEditor/` | Slicer 模块、编辑计算核心及合成演示 |
| `docs/` | 中英文算法及编辑器使用说明 |
| `tests/` | 数值测试和 Slicer 集成测试 |


