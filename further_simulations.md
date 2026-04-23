# Further Simulations Needed

## 1. Re-run `run_measure_extra.jl` with error bar support

The measurement script has been updated to compute per-k-point standard errors
(`hk_power_err`) for the structure factor. All existing data files need to be
regenerated to include this field.

### Commands to run

```bash
# For each (L, sigma) combination:
julia run_measure_extra.jl 64  0.5
julia run_measure_extra.jl 64  1.0
julia run_measure_extra.jl 64  2.0
julia run_measure_extra.jl 64  4.0
julia run_measure_extra.jl 256 0.5
julia run_measure_extra.jl 256 1.0
julia run_measure_extra.jl 256 2.0
julia run_measure_extra.jl 256 4.0
julia run_measure_extra.jl 512 0.5
julia run_measure_extra.jl 512 1.0
julia run_measure_extra.jl 512 2.0
julia run_measure_extra.jl 512 4.0
```

### Notes

- Default: 10000 thermalization + 100000 measurement sweeps
- After re-running, the data files will contain `hk_power_err` in addition to `hk_power_avg`
- The plotting script (`plot_histogram_chik.jl`) can then add error bars to the
  structure factor plot (Fig 4) once the new data is available
- The structure factor plot now only uses the x-direction (ky=0) instead of
  the full 2D radial average

## 2. Potential improvements

- Consider increasing measurement sweeps for better statistics on error bars
- Consider using blocking analysis for autocorrelation-corrected errors
