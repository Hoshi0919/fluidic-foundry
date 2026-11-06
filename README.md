# Fluidic Foundry (流控工场)

> **基于低雷诺数两相微流控水动力阻力切换与表面张力毛细效应的无晶体管流体逻辑与计算芯片**  
> *Microfluidic Bubble Logic, Hydrodynamic Resistance Inversion & Capillary Computing in Stokes Flow.*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.13](https://img.shields.io/badge/Python-3.13%2B-brightgreen.svg)](https://python.org)
[![Tests: 15 Passed](https://img.shields.io/badge/Tests-15%20Passed-emerald.svg)](tests/)
[![Live Demo](https://img.shields.io/badge/Demo-Interactive%20Canvas-cyan.svg)](https://hoshi0919.github.io/fluidic-foundry/)

---

## 1. 物理科学背景：破除低雷诺数线性对称性

在微流控芯片尺度（通道宽 $w \sim 50\,\mu\text{m}$，流速 $u \sim \text{mm/s}$），流体运动处于**极低雷诺数区（Stokes 蠕动流，$Re = \frac{\rho u D_h}{\mu} \ll 1$）**。

在单相牛顿流体中，纳维-斯托克斯方程退化为严格线性的 Stokes 方程：
$$0 = -\nabla P + \mu \nabla^2 \mathbf{u}$$

根据流体力学著名的**扇贝定理（Scallop Theorem，Purcell 1977）**，低雷诺数单相流具有严格的时间可逆性与线性叠加性。在没有外部机械电磁微阀门或活性部件的情况下，单相连续流体**无法产生非线性、双稳态或布尔逻辑运算**。

### 两相流引入非线性位元（Bits with Material Payload）

参考 Manu Prakash 与 Neil Gershenfeld 发表于 *Science* (2007, Vol. 315) 的奠基性成果，通过在油性连续相中引入不混溶的微水滴（或水性连续相中引入微气泡），系统的连续介质对称性被打破：

1. **信息载体物理化**：每个微液滴对应一个物理位元（Bit: $1$ 表示存在液滴，$0$ 表示纯连续相）；
2. **物质即信息**：与传统半导体中仅携带电子电荷不同，微液滴内部可包裹试剂、药物分子或细胞，**在完成逻辑寻路与控制运算的同时，直接执行微流体物料搬运与化学合成**；
3. **Bretherton 润滑膜非线性流阻反馈**：单个受限微栓进入管道时，由于弯月面高曲率拉普拉斯压力阶跃及微薄润滑膜的强烈剪切耗散，通道水动力阻力瞬间剧增：
   $$R_{\text{eff}} = R_0 + \kappa \cdot R_0 \cdot \frac{L_{\text{drop}}}{L_{\text{channel}}}$$
   这一流阻跃变反作用于全芯片压力场，使得后序液滴感知到下游拓扑流阻的变化，自发发生动力学转向。

---

## 2. 核心数学模型与逻辑门机理

### 2.1 矩形微通道泊肃叶流阻严格解析解
对于宽 $w$、高 $h$、长 $L$ 的矩形微截面通道：
$$R_0 = \frac{12 \mu L}{w h^3 \left[1 - \frac{192 h}{\pi^5 w} \sum_{n=1,3,...}^\infty \frac{1}{n^5} \tanh\left(\frac{n \pi w}{2 h}\right)\right]}$$

### 2.2 复合 AND / OR 非对称分流结 (Asymmetric Bifurcation Gate)
通过精密调节输出支路的几何截面，使得：
- **空载常态**：$R_{\text{OR}} < R_{\text{AND}} \implies Q_{\text{OR}} > Q_{\text{AND}}$（主流量偏向 OR 支路）；
- **首液滴进驻**：当首个液滴进入 OR 支路后，总阻力提升至 $R_{\text{OR}} + \Delta R > R_{\text{AND}}$；
- **流阻反转（Resistance Inversion）**：此时瞬时流量发生翻转：$Q_{\text{AND}} > Q_{\text{OR}}$；
- **布尔仲裁**：紧随其后的第二个液滴感知到 AND 支路成为最小阻力通道，被强迫分流进入 AND 支路！

| 输入 A | 输入 B | 输出 OR (支路 1) | 输出 AND (支路 2) | 水动力物理过程 |
| :---: | :---: | :---: | :---: | :--- |
| **0** | **0** | **0** | **0** | 无液滴，全场维持静态层流 |
| **1** | **0** | **1** | **0** | 单液滴顺应最小阻力，进入宽支路 OR |
| **0** | **1** | **1** | **0** | 单液滴顺应最小阻力，进入宽支路 OR |
| **1** | **1** | **1** | **1** | 首液滴占领 OR 支路引发流阻反转，次液滴被排斥挤入 AND 支路 |

### 2.3 NOT 反相器与侧向时钟偏转 (Clocked Inverter Gate)
连续注入时钟液滴流（Clock $CLK=1$）。当控制信号 $A=1$ 时，控制液滴在分流交叉口通过横向水动力跨阻桥施加压差偏转，将时钟液滴推挤入废液旁路（Dump），使得主输出端口在对应时钟周期内无液滴到达（输出 $\bar{A}=0$）。

### 2.4 1-Bit 全流体半加器 (Fluidic Half Adder)
同时实现：
- $\text{Sum} = A \oplus B$ (XOR 异或逻辑)
- $\text{Carry} = A \cdot B$ (AND 与逻辑)

---

## 3. 架构设计与代码组织

```
fluidic-foundry/
├── fluidic_foundry/
│   ├── __init__.py      # 模块统一入口与公共类型导出
│   ├── channel.py       # 矩形截面泊肃叶流阻解析解与 Bretherton 润滑膜液滴附加阻力模型
│   ├── network.py       # 基尔霍夫流体回路图求解器（高斯主元消元解全场压力与流速）
│   ├── simulator.py     # 拉格朗日液滴平流追踪、离散时间步推进与最大流速分流判据
│   └── gates.py         # 标准逻辑拓扑工厂（AND/OR 复合门、NOT 反相器、1-Bit 半加器）
├── tests/
│   └── test_fluidic_foundry.py  # 15 项物理强测试矩阵（流阻标度、连续性守恒、流阻反转、真值表闭合）
├── index.html           # 纯原生 HTML5 Canvas 交互式仿真器（双相流热力渐变、示踪粒子、Web Audio 音效）
└── README.md
```

---

## 4. 运行与验证

### 4.1 运行物理单元测试 (15 项全部通过)
```bash
PYTHONPATH=. uv run --with pytest pytest tests/ -v
```

输出：
```
============================= test session starts ==============================
collected 15 items

tests/test_fluidic_foundry.py::test_rectangular_channel_poiseuille_resistance_scaling PASSED [  6%]
tests/test_fluidic_foundry.py::test_droplet_excess_resistance_confinement PASSED [ 13%]
tests/test_fluidic_foundry.py::test_kirchhoff_flow_conservation_at_junction PASSED [ 20%]
tests/test_fluidic_foundry.py::test_hydrodynamic_ohms_law PASSED         [ 26%]
tests/test_fluidic_foundry.py::test_resistance_inversion_phase_transition PASSED [ 33%]
tests/test_fluidic_foundry.py::test_and_or_gate_truth_table_0_0 PASSED   [ 40%]
tests/test_fluidic_foundry.py::test_and_or_gate_truth_table_1_0 PASSED   [ 46%]
tests/test_fluidic_foundry.py::test_and_or_gate_truth_table_0_1 PASSED   [ 53%]
tests/test_fluidic_foundry.py::test_and_or_gate_truth_table_1_1 PASSED   [ 60%]
tests/test_fluidic_foundry.py::test_not_inverter_gate_uninhibited PASSED [ 66%]
tests/test_fluidic_foundry.py::test_half_adder_single_input_yields_sum PASSED [ 73%]
tests/test_fluidic_foundry.py::test_reynolds_and_capillary_regime_validity PASSED [ 80%]
tests/test_fluidic_foundry.py::test_half_adder_double_input_yields_carry PASSED [ 86%]
tests/test_fluidic_foundry.py::test_channel_geometry_validation PASSED   [ 93%]
tests/test_fluidic_foundry.py::test_network_singular_matrix_guard PASSED [100%]

============================== 15 passed in 0.10s ==============================
```

### 4.2 本地体验交互式芯片
在任意现代浏览器中打开 `index.html` 即可在无须任何构建工具的情况下体验纯原生高帧率微流控芯片交互。

---

## 5. 参考文献

1. **Manu Prakash, Neil Gershenfeld**, *Microfluidic Bubble Logic*, **Science**, Vol. 315, Issue 5813, pp. 832-835 (2007). DOI: [10.1126/science.1136907](https://doi.org/10.1126/science.1136907)
2. **E. M. Purcell**, *Life at low Reynolds number*, **American Journal of Physics**, 45, 3-11 (1977).
3. **F. P. Bretherton**, *The motion of long bubbles in tubes*, **Journal of Fluid Mechanics**, 10(2), 166-188 (1961).
4. **M. J. Fuerstman, P. Garstecki, G. M. Whitesides**, *Coding/Decoding and Reversibility of Droplet Trains in Microfluidic Networks*, **Science**, 315, 828-832 (2007).
5. **M. Schindler, A. Ajdari**, *Droplet traffic in microfluidic networks: A minimal model for network-level phenomena*, **Phys. Rev. Lett.**, 100, 044501 (2008).

---

## 作者与许可

- **作者**：Hoshi (`hoshi0919@outlook.com`)
- **协议**：MIT License
