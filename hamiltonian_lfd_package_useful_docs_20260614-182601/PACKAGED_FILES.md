# Hamiltonian / Hofstadter / LFD Package - Useful Docs Only

Created: 2026-06-14 18:26:01

## 保留的有用文档

- `README.md`: 项目基本入口。
- `TB_Hofstadter_LFD_mask_manual.md`: Hofstadter / DOS / LFD 工作流的主说明。
- `d3_two_band_peierls_landau_fan_manual_v5.md`: D3 two-band Peierls/Landau fan 的具体执行手册。
- `docs/model_versions/d3_slope_reversal_hamiltonian_v4.md`: 指定保留的 Hamiltonian v4 示例。

## 已移除的文档类型

- 重复 precheck 手册。
- 旧搜索过程手册。
- 过时或探索性质的 model notes。
- agent 本地记忆/环境说明。

## 核心代码

- `twse2_continuum_dos/peierls_hofstadter.py`: Hofstadter spectrum。
- `twse2_continuum_dos/dos.py`: DOS 和 density-axis LFD。
- `twse2_continuum_dos/d3_two_band.py`: D3 two-band Hamiltonian。
- `twse2_continuum_dos/d3_physical_reservoir.py`: 当前 physical reservoir Hamiltonian 与诊断。
- `examples/`: 可运行计算和绘图入口。
- `tests/`: 行为检测和最小示例。

## Packaged files

- `README.md`
- `pyproject.toml`
- `twse2_continuum_dos\__init__.py`
- `twse2_continuum_dos\peierls_hofstadter.py`
- `twse2_continuum_dos\tb_hofstadter_lfd_workflow.py`
- `twse2_continuum_dos\dos.py`
- `twse2_continuum_dos\d3_two_band.py`
- `twse2_continuum_dos\d3_physical_reservoir.py`
- `twse2_continuum_dos\d3_zero_field_diagnostics.py`
- `examples\run_tb_hofstadter_lfd_manual_quick.py`
- `examples\compute_d3_two_band_lfd_highres.py`
- `examples\render_d3_two_band_lfd_from_data.py`
- `examples\plot_d3_two_band_lfd.py`
- `examples\plot_d3_two_band_zero_field.py`
- `examples\plot_d3_physical_reservoir_diagnostics.py`
- `tests\test_peierls_hofstadter.py`
- `tests\test_tb_hofstadter_lfd_workflow.py`
- `tests\test_dos.py`
- `tests\test_d3_two_band.py`
- `tests\test_d3_physical_reservoir.py`
- `tests\test_d3_zero_field_diagnostics.py`
- `TB_Hofstadter_LFD_mask_manual.md`
- `d3_two_band_peierls_landau_fan_manual_v5.md`
- `docs\model_versions\d3_slope_reversal_hamiltonian_v4.md`
