# Upstream notes

These are things in OpenSWMM 6 (swmm6_rel @ 88329d6) that looked wrong or odd during the translation. The translation keeps every one as written, so results stay identical; this list is for the upstream team. Each entry quotes the code as it stands.

| # | Where | What | Can it change a run? |
|---|---|---|---|
| U1 | `math/MathExpr.cpp` `tokenize` | A lone `.` starts a "number", and `std::stod` throws `std::invalid_argument`. Nothing on the path from `ControlEngine::addExpression` catches it, so `EXPRESSION E1 = 2 + .` aborts the process (code-review deck IO-41_point). Code: `t.value = std::stod(s.substr(start, i - start));` | Yes: it aborts |
| U2 | `data/SubcatchData.hpp` `erase_at` | Does not erase `gage_name`, `snowpack_name`, `old_ponded_depth` or `washoff_vout1`; the list goes from `e(outlet_name);` to `e(gage);`. After an erase those four columns are out of step with the rest. | Only through the edit API (deleting a subcatchment) |
| U3 | `data/GageData.hpp` `erase_at` | Does not erase `station_id`, `rain_units`, `file_first_date`, `file_last_date`, `file_periods_precip` or `rain_series`; the header's own comment admits it. | Only through the edit API |
| U4 | `data/SubcatchData.hpp` | `sweep_now` and `sweep_last_date` are never sized, erased or shrunk by any method; `resize_coverage` sizes only `coverage` and `sweep_last_swept`. | Only if some later code indexes them unsized (to be checked when street sweeping is translated) |
| U5 | `data/SubcatchData.hpp` `resize` | `outlet_name`, `snowpack_name` and `gage_name` are `resize`d (old entries kept) while every other column is `assign`ed (cleared). | Only on a second `resize` of a populated store |
| U6 | `data/LinkData.hpp` `erase_at` | Does not erase `time_last_set`, `slot_volume`, `stat_slot_vol_dt`, `stat_vol_dt`, `stat_peak_slot_share` or `stat_time_slot_above`. | Only through the edit API |
| U7 | `data/LinkData.hpp` `reserve_to` | `r(stat_flow_class);` reserves `un`, where the array holds `un * N_FLOW_CLASSES` (compare `stat_lstep_bin_time`, reserved at `un * N_LSTEP_BINS`). | No (capacity only) |
| U8 | `data/LinkData.hpp` `resize` | `pump_curve_name.resize(un);` keeps old entries where the other columns are `assign`ed; `stat_pump_periods.assign(un, 0.0);` assigns a double to a `vector<long>`. | Only on a second `resize` of a populated store |
| U9 | `data/LinkSubtypes.hpp` `ConduitData::erase_at` | Leaves the five `gw_*` columns unerased (the doc comment says so). | Only through the edit API |
| U10 | `data/NodeSubtypes.hpp` `set_node_type` | Zeroes `storage_area` before the same-type early return, so setting STORAGE again on a storage node still clears its area: `if (t != NodeType::JUNCTION && ui < nodes.storage_area.size()) nodes.storage_area[ui] = 0.0;` (the differential shows it: area 7 becomes 0). | Only through the edit API |
| U11 | `data/NodeSubtypes.hpp` `set_node_type` | `shifted = true;` follows `switch (old)` even when `old` is JUNCTION with a row, a state that should not exist. | No |
| U12 | `data/InfraData.hpp` `InletUsageStore` | `add_row` does not grow the stat columns, but `erase_row` erases from them. | Only through the edit API |
| U13 | `data/NodeSubtypes.hpp` `erase_node` | `drop(storages) \|\| drop(outfalls) \|\| drop(dividers);` drops only the first matching row. | No, if a node has one subtype row, as it should |
