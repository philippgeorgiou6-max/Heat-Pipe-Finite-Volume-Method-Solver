# Heat-Pipe-Finite-Volume-Method-Solver
MATLAB finite-volume solver with GUI for steady 2D axisymmetric heat pipes. Solves compressible vapor flow, Brinkman–Darcy wick flow and wall conduction as one conjugate problem, coupled by evaporation and condensation. Independent implementation of Faghri &amp; Buchko (1991), validated against the paper.

<img width="800" height="398" alt="9dabdaa8-ee93-4d7a-97f5-731c2ffbe966" src="https://github.com/user-attachments/assets/73a7dcbe-44b6-4680-98f2-05c7616ebf86" />
<img width="800" height="444" alt="a62b3b99-b383-4882-bc33-a342398e8137" src="https://github.com/user-attachments/assets/378ae5c6-5fe0-408a-87d4-6ebffe52e9c0" />
<img width="800" height="414" alt="b3fd894e-78a2-4bab-a434-57fff7ef2174" src="https://github.com/user-attachments/assets/9e4dc5ab-9d51-4491-b8fc-fcd9c15ea09e" />
<img width="800" height="416" alt="a2919fe4-3572-4602-a66b-94c322675985" src="https://github.com/user-attachments/assets/b8e85f16-2410-445d-ad51-beb417b0b560" />
<img width="800" height="416" alt="a51c2c01-61bf-4ac9-8bf0-029d04e82c79" src="https://github.com/user-attachments/assets/cc113250-15cb-4907-a66c-6e56906a751b" />


# Heat Pipe CFD: conjugate 2D model of a copper–water heat pipe

![MATLAB](https://img.shields.io/badge/MATLAB-solver%20%2B%20GUI-orange)
![GNU Octave](https://img.shields.io/badge/GNU%20Octave-8.4%20tested-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-green)

A MATLAB finite-volume solver with a graphical user interface for the steady operation of a cylindrical heat pipe. It reproduces the numerical model of **A. Faghri and M. Buchko (1991)**, *Experimental and Numerical Analysis of Low-Temperature Heat Pipes With Multiple Heat Sources*, ASME J. Heat Transfer 113(3), 728–734, [doi:10.1115/1.2910624](https://doi.org/10.1115/1.2910624).

Vapor flow in the core, liquid flow through the wick and heat conduction in the wall are solved together as one conjugate problem. The three regions are coupled by evaporation and condensation at the liquid–vapor interface.

<p align="center">
  <img src="docs/images/gui_solver.png" alt="GUI after a run: live plots of temperature, velocity, axial temperatures, interface mass flux, pressures, convergence, power balance and thermal resistance" width="100%">
</p>
<p align="center"><em>The solver page of the GUI after a run of the paper's Case 1 (97 W).</em></p>

## Highlights

- **Full conjugate model:** a compressible, laminar vapor core, a Brinkman–Darcy liquid flow in the screen wick and conduction in the copper wall, all solved together.
- **Validated against the paper:** for Case 1 the computed wall temperature matches Faghri & Buchko's numerical solution within 0.2 K at the evaporator peak, along the transport section and along the condenser ([Validation](#validation)).
- **Six-step GUI:** geometry, boundary conditions, mesh, materials, numerics and solver. The live plots appear inside the GUI, and there are Stop and Save buttons.
- **Boundary conditions beyond the paper:** uniform or T(z) wall temperature, convection and q(z) at the evaporator and condenser, plus convection along the transport section. With these, the operating temperature is solved for.
- **Complete output:** every diagram as PDF and every data set as TXT, written to one results folder.
- **Scriptable:** `runHeatPipeSolver(...)` runs any case without the GUI.

## Contents

- [The physical model](#the-physical-model)
- [What the solution looks like](#what-the-solution-looks-like)
- [Validation](#validation)
- [Getting started](#getting-started)
- [The GUI](#the-gui)
- [Boundary conditions](#boundary-conditions)
- [Output](#output)
- [Repository structure](#repository-structure)
- [Limitations](#limitations)
- [References](#references)
- [License](#license)

## The physical model

### The heat pipe (paper Table 1)

| Parameter | Value |
|---|---|
| Length | 1000 mm |
| Wall | copper, outer diameter 25.4 mm, inner diameter 22.0 mm, k = 339 W/(m K) |
| Wick | two wraps of 50-mesh copper screen (N = 1970 m⁻¹, wire diameter 0.178 mm), 0.712 mm thick |
| Wick properties (Chi 1976) | porosity 0.711, permeability 1.12 × 10⁻⁹ m², capillary radius 0.254 mm, k_eff = 1.20 W/(m K) |
| Working fluid | water; operating temperature T_o = 75.7 °C (Case 1) |
| Evaporator | 63.5 mm, from z = 20 mm to 83.5 mm, 97 W |
| Transport section | 596.5 mm, adiabatic |
| Condenser | 300 mm, from z = 680 mm to 980 mm; removes the input as a uniform heat flux |
| End sections | 20 mm at each end, adiabatic |

The paper's pipe has four evaporators. This project models one, which matches the paper's Case 1, where only evaporator 1 is heated.

### Governing equations

The model is steady, laminar and axisymmetric (coordinates r, z; velocity **V** = (v, w)).

**Vapor core** (Eqs. 1–4): compressible Navier–Stokes and energy equations with pressure work and viscous dissipation; the vapor is an ideal gas.

$$\nabla\cdot(\rho_v\mathbf{V})=0,\qquad \rho_v(\mathbf{V}\cdot\nabla)\mathbf{V}=-\nabla p+\mu_v\left[\nabla^2\mathbf{V}+\tfrac{1}{3}\nabla(\nabla\cdot\mathbf{V})\right]$$

$$\rho_v c_{p,v}\,\mathbf{V}\cdot\nabla T=k_v\nabla^2 T+\mathbf{V}\cdot\nabla p+\mu_v\Phi,\qquad \rho_v=\frac{p}{RT}$$

**Wick** (Eqs. 5–9): the liquid-saturated screen is a porous medium (porosity ε, permeability K, intrinsic velocity **V**) with Brinkman–Darcy momentum including inertia. The effective conductivity k_eff follows Chi (1976), paper Eq. 9.

$$\nabla\cdot(\rho_l\varepsilon\mathbf{V})=0,\qquad \rho_l(\mathbf{V}\cdot\nabla)\mathbf{V}=-\nabla p+\mu_l\nabla^2\mathbf{V}-\frac{\mu_l\varepsilon}{K}\mathbf{V},\qquad \rho_l c_{p,l}\,\varepsilon\mathbf{V}\cdot\nabla T=k_\mathrm{eff}\nabla^2 T$$

**Wall** (Eq. 10): heat conduction, $\nabla^2 T=0$.

**Liquid–vapor interface** r = a: the temperature is the saturation temperature of the local vapor pressure (Eq. 19). The heat conducted to the interface evaporates or condenses the liquid (Eq. 16), and there is no slip (Eq. 13).

$$T_i=\left[\frac{1}{T_o}-\frac{R}{h_{fg}}\ln\frac{p_v}{p_o}\right]^{-1},\qquad \dot m''=\frac{1}{h_{fg}}\left(k_\mathrm{eff}\frac{\partial T_l}{\partial r}-k_v\frac{\partial T_v}{\partial r}\right),\qquad v_v=-\frac{\dot m''}{\rho_v},\qquad v_l=-\frac{\dot m''}{\rho_l\varepsilon},\qquad w=0$$

**Outer wall** r = c (Eq. 20): q = Q/A_e into the evaporator, q = −Q/A_c out of the condenser and q = 0 elsewhere. The end caps are adiabatic with no slip, and the axis is a symmetry line.

Water properties come from IAPWS-based correlations, evaluated per region at the operating temperature T_o. The saturation pressure p_o = p_sat(T_o) and the local latent heat h_fg(T_i) follow the same correlations.

### Numerical method

- Staggered finite-volume grid in r and z, with hybrid differencing by default; power law and upwind are also available.
- SIMPLEC pressure–velocity coupling in the vapor core and in the wick; SIMPLE is also available. The paper used Spalding's (1980) finite-domain method with the SIMPLEST practice.
- One energy equation for wall, wick and vapor, with harmonic-mean conductivities at the region boundaries. The interface is an internal face held at T_i from Eq. 19, and the evaporation rate follows from the jump of the conductive heat flux across it (Eq. 16).
- Each outer iteration consists of the flow iterations in vapor and wick, then the energy equation with the interface conditions, then the update of the interface mass flux. Iteration stops when all convergence criteria (default 10⁻⁶) hold in the same outer iteration, and the residuals are checked again on the final fields.
- The default grid has 335 axial × (24 vapor + 16 wick + 8 wall) = 16,080 cells. Cells are clustered at both faces of the wick to resolve its Brinkman boundary layers. Case 1 converges in 84 outer iterations, about a minute.

## What the solution looks like

<p align="center">
  <img src="docs/images/fields_case1.png" alt="Temperature field in wall, wick and vapor core, and vapor streamlines from evaporator to condenser" width="100%">
</p>

For Case 1 (97 W):

- **Wall and wick:** the evaporator heats the wall to 86.8 °C, and the condenser cools it to 72.8 °C.
- **Vapor core:** the temperature stays at T_o = 75.7 °C. Along the interface it varies by only 0.2 mK, because the vapor pressure drop is just 0.36 Pa.
- **Vapor flow:** vapor leaves the heated wick and reaches up to 1.0 m/s on the axis (axial Reynolds number below 230, so the flow is laminar). It condenses along the condenser.
- **Liquid return:** the liquid flows back through the wick at up to 1.4 mm/s.

The radial direction is stretched in both plots, so that each region is visible.

## Validation

### Comparison with Faghri & Buchko (1991), Case 1

<p align="center">
  <img src="docs/images/validation_case1.png" alt="Wall and vapor temperature along the pipe: this model, the paper's model and the experiment" width="100%">
</p>

The paper's computed curves and the measurements were digitized from its Fig. 4, with an estimated accuracy of ±0.2 K. The data are in [`docs/validation/faghri_buchko_1991_fig4_digitized.csv`](docs/validation/faghri_buchko_1991_fig4_digitized.csv).

| | This model | Paper's model | Experiment |
|---|---|---|---|
| Wall, evaporator peak | 86.8 °C | 86.7 °C | 85.6 °C |
| Wall peak − vapor | 11.1 K | 11.1 K | 10.0 K |
| Wall, condenser | 72.8 °C | 72.7 °C | 74.2–75.1 °C |
| Vapor | 75.70 °C | 75.6 °C | 75.4–75.7 °C |

- **Against the paper's model:** this model reproduces the paper's computed wall temperature within 0.2 K at the evaporator peak, along the transport section and along the condenser (median difference 0.08 K). On the steep flanks of the peak the differences reach 0.4 K, which corresponds to an axial shift of about 2 mm, within the accuracy of the digitized curve.
- **Against the experiment:** both models predict an evaporator peak about 1 K too high and a condenser wall 1.5–2 K too low. The paper attributes the condenser difference to the uniform heat flux assumed at the condenser (Eq. 20), which does not match the water calorimeter used in the test, and to the cooling water affecting the condenser thermocouples.
- **Wick thickness:** Table 1 lists a 0.712 mm wick and a 20.5 mm vapor core, but those diameters imply 0.75 mm. With 0.712 mm (the default) the model matches the paper's curve. With 0.75 mm the evaporator peak comes out 0.5 K higher, shown dotted in the right panel above. The 0.75 mm variant is selected with `wickThicknessSource = 'diameters'`.

### Grid independence

<p align="center">
  <img src="docs/images/grid_study.png" alt="Wall temperature on four grids in the evaporator and condenser regions" width="100%">
</p>

| Grid | Axial × (vapor + wick + wall) | Cells | Outer iterations | Wall peak − vapor | Thermal resistance | Wick pressure drop |
|---|---|---:|---:|---:|---:|---:|
| paper50 | 50 × (10 + 5 + 5) | 1,000 | 19 | 11.21 K | 0.134 K/W | 248 Pa |
| paper100 (the paper's grid) | 100 × (20 + 5 + 10) | 3,500 | 55 | 11.08 K | 0.131 K/W | 248 Pa |
| default | 335 × (24 + 16 + 8) | 16,080 | 84 | 11.13 K | 0.130 K/W | 262 Pa |
| fine | 669 × (48 + 32 + 16) | 64,224 | 249 | 11.14 K | 0.130 K/W | 263 Pa |

- **Temperatures:** they agree within 0.1 K on all four grids. The default grid differs from the fine grid by less than 0.1 % in the temperatures and the thermal resistance, and by 0.4 % in the pressure drops.
- **Wick pressure drop on the paper's grids:** it is 6 % lower there. Five uniform wick cells cannot resolve the 0.04 mm Brinkman layer at the wick faces; the clustered default grid can.
- **Thermal resistance:** defined as (mean evaporator wall − mean condenser wall temperature) / Q.

### Checks against one-dimensional theory

| Quantity | 1-D estimate | This model |
|---|---|---|
| Wall − vapor: peak / evaporator end cap / condenser | 11.11 / 4.72 / −2.71 K (copper wall as a fin on the wick) | 11.14 / 4.72 / −2.71 K |
| Wick pressure drop | 264 Pa (Darcy–Brinkman; Darcy alone 235 Pa) | 262 Pa |
| Vapor pressure drop | 0.342 Pa (Poiseuille) | 0.363 Pa (+6 % from inertia, radial Reynolds number 7.6) |
| Variation of the interface temperature | 0.221 mK (Clausius–Clapeyron) | 0.221 mK |
| Evaporation rate | Q / h_fg = 4.183 × 10⁻⁵ kg/s | 4.183 × 10⁻⁵ kg/s |

Mass and energy are conserved to the convergence tolerance: heat in equals heat out, and everything that evaporates condenses.

## Getting started

**Requirements:** MATLAB (written for R2016b or later; no toolboxes). Solver and GUI were also tested in GNU Octave 8.4.

```matlab
RUN_ME            % opens the GUI
```

Without the GUI:

```matlab
results = runHeatPipeSolver();                                              % paper Case 1
results = runHeatPipeSolver(struct('evaporatorPowers_W',150,'T_o_C',75));
results = runHeatPipeSolver(struct('condenserBCType','convection','condenserH_W_m2K',2000));
results = runHeatPipeSolver([], 'paper100');                                % grid preset
results = runHeatPipeSolver([], 'default', struct('alpha_u',0.5,'tolMass',1e-5));  % numerics
```

All case inputs are listed in [`src/input/inputCleanHeatPipeCase.m`](src/input/inputCleanHeatPipeCase.m), and the numerical settings and convergence criteria in [`src/input/inputCleanHeatPipeNumerics.m`](src/input/inputCleanHeatPipeNumerics.m). A script run writes the same results folder as the GUI. When a display is available, it also shows live plots and every diagram in its own window.

## The GUI

The GUI leads through six steps, and "Next step" checks each step before it moves on.

<table>
  <tr>
    <td width="33%"><img src="docs/images/gui_geometry.png" alt="Step 1: geometry"></td>
    <td width="33%"><img src="docs/images/gui_boundary_conditions.png" alt="Step 2: boundary conditions"></td>
    <td width="33%"><img src="docs/images/gui_mesh.png" alt="Step 3: mesh"></td>
  </tr>
  <tr>
    <td align="center">1 Geometry</td>
    <td align="center">2 Boundary conditions</td>
    <td align="center">3 Mesh</td>
  </tr>
  <tr>
    <td width="33%"><img src="docs/images/gui_materials.png" alt="Step 4: materials"></td>
    <td width="33%"><img src="docs/images/gui_numerics.png" alt="Step 5: numerics"></td>
    <td width="33%"><img src="docs/images/gui_solver.png" alt="Step 6: solver and results"></td>
  </tr>
  <tr>
    <td align="center">4 Materials</td>
    <td align="center">5 Numerics</td>
    <td align="center">6 Solver &amp; results</td>
  </tr>
</table>

| Step | Content |
|---|---|
| 1 Geometry | Lengths of the unheated sections, the evaporator and the condenser, plus the three diameters. "Define geometry" draws the pipe with its regions and dimensions. The transport section takes the remaining length. |
| 2 Boundary conditions | Choose a type for the evaporator, the transport section and the condenser. Only the fields of the selected type are shown, and a preview plots what the solver will impose. |
| 3 Mesh | Radial cells and spacing per region, axial cell size, and presets (including the paper's grids). "Start meshing" draws the mesh in full and as a detail. |
| 4 Materials | Vapor, wick (liquid plus porous structure; k_eff from Eq. 9 or user-defined) and wall (material presets or custom values). By default the water properties follow the operating temperature. |
| 5 Numerics | Pressure–velocity coupling, convection scheme, wick momentum form, start-up, inner iterations and under-relaxation factors. Also the convergence criteria: maximum outer iterations and tolerances for mass, momentum, temperature change, interface mass-flux change, interface closure, energy residual and heat balance. |
| 6 Solver & results | Initial temperature and live-plot interval, with Start simulation, Stop and Save results. |

The geometry and mesh views stretch the radial direction: the vapor core, wick and wall take 50, 25 and 25 % of the height, and the tick labels show the true radii.

## Boundary conditions

| Zone | Types |
|---|---|
| Evaporator | input power (paper), uniform temperature, T(z), convection, q(z) |
| Transport section | adiabatic (paper), convection |
| Condenser | remove the input power (paper Eq. 20), power removed, uniform temperature, T(z), convection, q(z) |

- **Function inputs:** T(z) and q(z) are MATLAB expressions in `z` (axial position, m), `s` (distance from the start of the zone, m) and `L` (zone length, m), for example `85 + 3*sin(pi*s/L)`. Temperatures are in °C and heat fluxes in W/m².
- **Temperature level:** if only heat fluxes are prescribed, the operating temperature T_o sets the level, as in the paper. With a temperature or convection condition, T_o is solved for so that everything that evaporates also condenses. The water properties then follow T_o during the run.
- **Validation status:** the paper's flux conditions are the validated configuration. The other types are extensions. They are checked for convergence and energy balance, but not validated against measurements.

## Output

"Save results" writes `./heatpipe_results`, after clearing the folder, so it always holds exactly one simulation.

| File | Content |
|---|---|
| `*.pdf` | Geometry, mesh, boundary conditions, temperature and velocity fields, axial temperatures, interface mass flux, pressures, convergence, power in/out, thermal resistance, and an overview (`diagnostics.pdf`) |
| `*.txt` | The data of every plot (tab-separated, `#` header lines), plus `summary.txt`, `geometry.txt`, `mesh.txt`, `boundary_conditions.txt` and `materials.txt` |
| `heatpipe_results.mat` | The complete results structure |

## Repository structure

| Path | Content |
|---|---|
| `RUN_ME.m` | Starts the GUI |
| `HeatPipeGUI.m` | The GUI |
| `src/runHeatPipeSolver.m` | Script entry point: solve, post-process, export, plot |
| `src/runPaperModelSolver.m` | The coupled outer iteration |
| `src/input` | Case, materials and numerical settings |
| `src/mesh` | Block-structured axisymmetric mesh, regions and zones |
| `src/flow`, `src/physics`, `src/numerics` | Momentum, pressure correction, energy, interface coupling, boundary conditions and water properties |
| `src/output` | Derived quantities, thermal resistance and export |
| `src/plotting` | Diagrams, geometry/mesh/BC drawings and the live monitor |
| `docs` | README images and the digitized validation data |

## Limitations

The model shares the assumptions of the paper:

- **Flow and state:** steady, laminar and axisymmetric, with no gravity. The vapor is an ideal gas, and the interface is in thermodynamic equilibrium (no kinetic resistance).
- **Liquid:** the wick is always saturated. The meniscus and the capillary pressure are not coupled to the flow, so dry-out and the capillary limit are not predicted; the capillary check in the output is a post-processed estimate.
- **Properties:** constant within each region, at the operating temperature.
- **Condenser:** in the paper configuration it is a uniform heat flux, which is the main reason for the difference from the measured condenser temperatures.
- **Scope:** only one evaporator is set up. The paper's multi-evaporator cases need the four-heater layout.

## References

- A. Faghri, M. Buchko (1991). Experimental and Numerical Analysis of Low-Temperature Heat Pipes With Multiple Heat Sources. *ASME Journal of Heat Transfer* 113(3), 728–734. [doi:10.1115/1.2910624](https://doi.org/10.1115/1.2910624)
- S. W. Chi (1976). *Heat Pipe Theory and Practice: A Sourcebook*. Hemisphere, Washington.
- S. V. Patankar (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere, Washington.
- J. P. Van Doormaal, G. D. Raithby (1984). Enhancements of the SIMPLE method for predicting incompressible fluid flows. *Numerical Heat Transfer* 7, 147–163.

## License

MIT License, see [`LICENSE`](LICENSE). The physical model is that of Faghri & Buchko (1991); please cite the paper when you use the model.
