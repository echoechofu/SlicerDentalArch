"""Interactive dental width and arch revision editor for 3D Slicer. BSD-3-Clause."""
import json
import logging
from pathlib import Path
import qt
import slicer
from slicer.ScriptedLoadableModule import ScriptedLoadableModule, ScriptedLoadableModuleWidget, ScriptedLoadableModuleLogic, ScriptedLoadableModuleTest
from DentalArchEditorLib.session import EditorSession
from DentalArchEditorLib import demo, geometry, stateio


class DentalArchEditor(ScriptedLoadableModule):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent.title = "Dental Arch Editor / 牙弓修订"
        self.parent.categories = ["Dentistry"]
        self.parent.contributors = ["DentalArchEditor contributors"]
        self.parent.helpText = "Edit existing crown-width landmarks and planar B-spline reference arches. Save new revisions without writing the source files. Native coordinate units; not a clinically validated diagnostic system."
        self.parent.acknowledgementText = "Built with 3D Slicer, NumPy, SciPy and VTK. BSD-3-Clause extension."


class DentalArchEditorLogic(EditorSession, ScriptedLoadableModuleLogic):
    pass


class DentalArchEditorWidget(ScriptedLoadableModuleWidget):
    def setup(self):
        super().setup()
        self.parent.setStyleSheet("QLabel, QComboBox, QPushButton, QCheckBox, QTableWidget { font-family: sans-serif; font-size: 12px; }")
        self._refreshing = False
        self.logic = None
        title = qt.QLabel("<b>牙弓与拥挤度修订 / Dental Arch Editor</b>")
        title.setWordWrap(True)
        self.layout.addWidget(title)
        intro = qt.QLabel("打开 DentalArchScene.mrb 或修订场景即可继续编辑。旧版标注场景需导入对应 JSON。\nOpen a generated scene or revision to edit. Legacy scenes require baseline JSON import.")
        intro.setWordWrap(True)
        self.layout.addWidget(intro)
        self.importButton = qt.QPushButton("导入原算法数据 / Import baseline")
        self.demoButton = qt.QPushButton("加载合成示例 / Synthetic demo")
        self.layout.addWidget(self.importButton)
        self.layout.addWidget(self.demoButton)
        self.toothBox = qt.QComboBox()
        form = qt.QFormLayout()
        form.addRow("当前牙齿 / Tooth", self.toothBox)
        self.crownBox = slicer.qMRMLNodeComboBox()
        self.crownBox.nodeTypes = ["vtkMRMLModelNode"]
        self.crownBox.noneEnabled = True
        self.crownBox.addEnabled = False
        self.crownBox.removeEnabled = False
        self.crownBox.setMRMLScene(slicer.mrmlScene)
        form.addRow("所属牙冠 / Crown", self.crownBox)
        self.layout.addLayout(form)
        self.modeBox = qt.QComboBox()
        self.modeBox.addItem("调整接触测量点 / Edit M-D points", "points")
        self.modeBox.addItem("调整牙弓手柄 / Edit arch handles", "arch")
        self.modeBox.setToolTip("黄色手柄控制弓形，通常不在弓线上；两端绑定第一磨牙近中点。\nYellow handles control shape and may lie off the curve. Endpoints follow first-molar mesial points.")
        self.layout.addWidget(self.modeBox)
        self.snapBox = qt.QCheckBox("吸附到所属牙冠 / Snap to own crown")
        self.snapBox.setChecked(True)
        self.compareBox = qt.QCheckBox("显示原算法对照 / Show baseline")
        self.layout.addWidget(self.snapBox)
        self.layout.addWidget(self.compareBox)
        self.totalLabel = qt.QLabel("尚未载入 / No session")
        self.totalLabel.setStyleSheet("font-size: 20px; font-weight: bold; padding: 8px;")
        self.totalLabel.setWordWrap(True)
        self.layout.addWidget(self.totalLabel)
        self.metrics = qt.QTableWidget(4, 4)
        self.metrics.setHorizontalHeaderLabels(["测量", "原始", "当前", "变化 Δ"])
        self.metrics.verticalHeader().hide()
        self.metrics.setEditTriggers(qt.QAbstractItemView.NoEditTriggers)
        self.metrics.setSelectionMode(qt.QAbstractItemView.NoSelection)
        self.metrics.setMaximumHeight(170)
        self.layout.addWidget(self.metrics)
        self.widthLabel = qt.QLabel("")
        self.widthLabel.setWordWrap(True)
        self.layout.addWidget(self.widthLabel)
        self.statusLabel = qt.QLabel("")
        self.statusLabel.setWordWrap(True)
        self.layout.addWidget(self.statusLabel)
        row = qt.QHBoxLayout()
        self.undoButton = qt.QPushButton("撤销 / Undo")
        self.redoButton = qt.QPushButton("重做 / Redo")
        row.addWidget(self.undoButton)
        row.addWidget(self.redoButton)
        self.layout.addLayout(row)
        self.resetButton = qt.QPushButton("恢复当前牙原始点位 / Reset tooth")
        self.refitButton = qt.QPushButton("按当前点重新拟合牙弓 / Refit arch")
        self.saveButton = qt.QPushButton("保存为新修订版 / Save new revision")
        for button in [self.resetButton, self.refitButton, self.saveButton]:
            self.layout.addWidget(button)
        note = qt.QLabel("手动修订模式：改内部牙点只更新牙宽；第一磨牙边界会联动。重拟合按钮将替换当前弓形，可撤销。\nNative units. C = R − A: positive shortage, negative surplus. Synthetic demo is not clinical data.")
        note.setWordWrap(True)
        self.layout.addWidget(note)
        self.layout.addStretch(1)
        self.logic = DentalArchEditorLogic(self.refresh)
        self.importButton.connect("clicked()", self.onImport)
        self.demoButton.connect("clicked()", lambda: self.runAction(lambda: demo.load_demo(self.logic)))
        self.toothBox.connect("currentIndexChanged(int)", self.onSelection)
        self.modeBox.connect("currentIndexChanged(int)", self.onSelection)
        self.crownBox.connect("currentNodeChanged(vtkMRMLNode*)", self.onCrown)
        self.snapBox.connect("toggled(bool)", self.onSnap)
        self.compareBox.connect("toggled(bool)", self.onCompare)
        self.undoButton.connect("clicked()", lambda: self.runAction(self.logic.undo))
        self.redoButton.connect("clicked()", lambda: self.runAction(self.logic.redo))
        self.resetButton.connect("clicked()", lambda: self.runAction(lambda: self.logic.reset_tooth(self.toothBox.currentText)))
        self.refitButton.connect("clicked()", self.onRefit)
        self.saveButton.connect("clicked()", self.onSave)
        self.refresh(self.logic.result)

    def cleanup(self):
        if self.logic:
            self.logic.cleanup()
        super().cleanup()

    def enter(self):
        if self.logic:
            self.logic.restore()
            self.refresh(self.logic.result)

    def runAction(self, action):
        try:
            return action()
        except Exception as error:
            logging.exception("DentalArchEditor action failed")
            slicer.util.errorDisplay(str(error))

    def onImport(self):
        md = qt.QFileDialog.getOpenFileName(self.parent, "近远中测量 JSON / MD measurements", "", "JSON (*.json)")
        if not md:
            return
        arch = qt.QFileDialog.getOpenFileName(self.parent, "参考牙弓 JSON / Reference arch", str(Path(md).parent), "JSON (*.json)")
        if arch:
            self.runAction(lambda: self.logic.import_baseline(md, arch))

    def onSelection(self, *_):
        if self._refreshing or not self.logic or not self.logic.parameter_node:
            return
        mode = "arch" if self.modeBox.currentIndex == 1 else "points"
        self.logic.set_edit_mode(mode, self.toothBox.currentText)
        self.refresh(self.logic.result)

    def onCrown(self, node):
        if not self._refreshing and self.logic and self.logic.parameter_node:
            self.runAction(lambda: self.logic.set_crown(self.toothBox.currentText, node))

    def onSnap(self, checked):
        if not self._refreshing and self.logic and self.logic.parameter_node:
            self.logic.parameter_node.SetParameter("SnapEnabled", "1" if checked else "0")

    def onCompare(self, checked):
        if not self._refreshing and self.logic and self.logic.parameter_node:
            self.logic.show_baseline(checked)

    def onRefit(self):
        self.statusLabel.setText("正在重新拟合 / Recomputing…")
        self.refitButton.setEnabled(False)
        slicer.app.processEvents()
        self.runAction(self.logic.refit)
        self.refresh(self.logic.result)

    def onSave(self):
        folder = qt.QFileDialog.getExistingDirectory(self.parent, "选择修订版保存目录 / Revision parent directory")
        if folder:
            saved = self.runAction(lambda: self.logic.save_revision(folder))
            if saved:
                self.statusLabel.setText("已另存 / Saved new revision:\n" + str(saved))

    def refresh(self, result):
        if self._refreshing or self.logic is None:
            return
        self._refreshing = True
        try:
            active = self.logic.parameter_node is not None and self.logic.state is not None
            for widget in [self.toothBox, self.crownBox, self.modeBox, self.snapBox, self.compareBox, self.resetButton, self.undoButton, self.redoButton]:
                widget.setEnabled(active)
            self.importButton.setEnabled(not active)
            self.demoButton.setEnabled(not active)
            self.saveButton.setEnabled(active and bool(result.get("valid")))
            self.refitButton.setEnabled(active and bool(result.get("valid")))
            if active:
                names = list(self.logic.state["teeth"])
                if [self.toothBox.itemText(i) for i in range(self.toothBox.count)] != names:
                    self.toothBox.clear()
                    self.toothBox.addItems(names)
                selected = self.logic.parameter_node.GetParameter("SelectedTooth") or names[0]
                self.toothBox.setCurrentIndex(max(0, self.toothBox.findText(selected)))
                self.crownBox.setCurrentNode(self.logic.parameter_node.GetNodeReference("Crown_" + selected))
                self.modeBox.setCurrentIndex(1 if self.logic.parameter_node.GetParameter("EditMode") == "arch" else 0)
                self.snapBox.setChecked(self.logic.parameter_node.GetParameter("SnapEnabled") == "1")
                baseline = geometry.evaluate(self.logic.baseline)
            else:
                baseline = {}
            labels = [("牙宽总和 R", "R"), ("参考弧长 A", "A"), ("拥挤量 C=R−A", "C"), ("空间余额 D=A−R", "D")]
            for row, (label, key) in enumerate(labels):
                delta = result[key] - baseline[key] if result.get(key) is not None and baseline.get(key) is not None else None
                values = [label, "—" if baseline.get(key) is None else f"{baseline[key]:.3f}", "—" if result.get(key) is None else f"{result[key]:.3f}", "—" if delta is None else f"{delta:+.3f}"]
                for col, value in enumerate(values):
                    self.metrics.setItem(row, col, qt.QTableWidgetItem(value))
            self.metrics.resizeColumnsToContents()
            if result.get("valid"):
                c = result["C"]
                self.totalLabel.setText((f"拥挤不足 / Shortage  {c:.3f}" if c > 0 else f"净余隙 / Surplus  {-c:.3f}") + "\n原始坐标单位 / native units")
                selected = self.toothBox.currentText
                if selected in self.logic.state["included_teeth"]:
                    scope = "计入牙宽总和 R / Included in R"
                elif selected in [b["tooth"] for b in self.logic.state["boundaries"]]:
                    scope = "宽度不计入 R；近中点控制弓端点 / Boundary tooth"
                else:
                    scope = "宽度不计入本次 R / Outside measurement range"
                self.widthLabel.setText(f"{selected}  宽度 / Width: {result['widths'][selected]:.3f}\n{scope}")
            else:
                self.totalLabel.setText("待修正 / Measurement unavailable")
                self.widthLabel.setText("")
            messages = list(result.get("errors", [])) + list(result.get("warnings", []))
            self.statusLabel.setText("\n".join(messages[:5]) + (f"\n… {len(messages)} 条提示 / notices" if len(messages) > 5 else "") if messages else "实时结果已更新 / Live values updated")
            self.undoButton.setEnabled(active and (self.logic._cursor > 0 or self.logic._structural_invalid))
            self.redoButton.setEnabled(active and self.logic._cursor + 1 < len(self.logic._history))
        finally:
            self._refreshing = False


class DentalArchEditorTest(ScriptedLoadableModuleTest):
    def runTest(self):
        slicer.mrmlScene.Clear(0)
        logic = DentalArchEditorLogic()
        try:
            result = demo.load_demo(logic)
            self.assertTrue(result["valid"])
            before = result["A"]
            controls = logic.state["arch"]["control_points_xy"]
            logic.set_arch_control(3, [controls[3][0], controls[3][1] + 1])
            self.assertGreater(logic.result["A"], before)
            logic.undo()
            self.assertAlmostEqual(logic.result["A"], before, places=8)
        finally:
            logic.cleanup()
