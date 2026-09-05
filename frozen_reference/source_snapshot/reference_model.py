"""Conditional ALBA reference plant, physical states x, v, V, Q (SI units).

Adds liquid inertia and a smooth-pipe Reynolds-dependent friction law to the
unchanged Stage-4 mechanics/actuator. No waterway or material refitting here.
The algebraic pressure is not a new independent bladder compliance state.
No voltage-duty DC motor, collapse law, post-slack or reverse/air flow model.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "baseline_step4"))
import step4_actuator as step4

Command = step4.Command
load_plant = step4.load_plant
load_actuator = step4.load_actuator


def load_reference(p, **overrides):
    h = json.loads((ROOT / "reference_config.json").read_text(encoding="utf-8"))
    unknown = set(overrides) - set(h)
    if unknown:
        raise ValueError("Unknown reference parameter: " + ", ".join(sorted(unknown)))
    h.update(overrides)
    area = np.pi * p["pipe_diameter"]**2 / 4
    length = sum(p["pipe_segments"])
    h["pipe_area_m2"] = area
    h["pipe_length_m"] = length
    h["R_pipe_laminar_Pa_s_m3"] = 128*p["water_viscosity"]*length/(np.pi*p["pipe_diameter"]**4)
    h["R_residual_Pa_s_m3"] = p["Rh"] - h["R_pipe_laminar_Pa_s_m3"]
    h["Ih_Pa_s2_m3"] = h["hydraulic_inertia_scale"]*p["rho"]*length/area
    h["linear_hydraulic_time_constant_s"] = h["Ih_Pa_s2_m3"] / p["Rh"]
    if not (h["Ih_Pa_s2_m3"] > 0 and h["R_residual_Pa_s_m3"] >= 0):
        raise ValueError("Positive inertia and nonnegative residual resistance are required.")
    if not (h["Re_turbulent"] > h["Re_laminar"] > 1 and h["pipe_absolute_roughness_m"] >= 0):
        raise ValueError("Invalid transition thresholds or roughness.")
    if h["additional_component_K_Pa_s2_m6"] < 0:
        raise ValueError("A passive hydraulic loss coefficient cannot be negative.")
    if h["additional_component_K_Pa_s2_m6"] and not h["additional_component_K_authorized"]:
        raise ValueError("Nonzero component K needs explicit calibrated/scenario authorization.")
    if h["material_relaxation_enabled"]:
        raise ValueError("Material relaxation is not identified or implemented in this version.")
    return h


def hydraulic_loss(Q, p, h):
    """Odd, passive pressure-loss curve; exactly Rh*Q in the laminar limit.

    Do NOT add Rh*Q and the full pipe loss: Rh already includes the pipe.
    R_residual is a conditional decomposition of the fitted total, not an
    independently measured resistance of the valve/flowmeter.
    """
    q = np.asarray(Q, dtype=float)
    Re = p["rho"]*np.abs(q)*p["pipe_diameter"]/(p["water_viscosity"]*h["pipe_area_m2"])
    laminar = h["R_pipe_laminar_Pa_s_m3"] * q
    weight = np.clip((Re-h["Re_laminar"])/(h["Re_turbulent"]-h["Re_laminar"]), 0., 1.)
    weight = weight**2*(3-2*weight)
    # Never evaluate turbulent logarithms at zero or in their invalid low-Re region.
    re_safe = np.maximum(Re, h["Re_laminar"])
    log_arg = 6.9/re_safe + (h["pipe_absolute_roughness_m"]/(3.7*p["pipe_diameter"]))**1.11
    f_turb = (-1.8*np.log10(log_arg))**-2
    turbulent = f_turb*h["pipe_length_m"]/p["pipe_diameter"]*p["rho"]/(2*h["pipe_area_m2"]**2)*q*np.abs(q)
    pipe = laminar + weight*(turbulent-laminar) if h["nonlinear_pipe"] else laminar
    residual = h["R_residual_Pa_s_m3"]*q
    optional = h["additional_component_K_Pa_s2_m6"]*q*np.abs(q)
    return dict(dp_loss=pipe+residual+optional, dp_pipe=pipe, dp_residual=residual,
                dp_optional=optional, Re=Re, transition_weight=weight)


def steady_flow(drive_Pa, p, h):
    """Forward-flow static limit, for diagnostics, not a prescribed-flow input."""
    if drive_Pa <= 0:
        return 0.
    upper = max(drive_Pa/p["Rh"], 1e-9)
    while hydraulic_loss(upper, p, h)["dp_loss"] < drive_Pa:
        upper *= 2
    return brentq(lambda q: hydraulic_loss(q, p, h)["dp_loss"]-drive_Pa, 0., upper,
                  xtol=1e-16, rtol=1e-12)


def hydraulic_slope(Q, p, h):
    """Analytic d(delta-p)/dQ; used to integrate fast liquid inertia stably."""
    qabs = abs(float(Q)); Rlam = h["R_pipe_laminar_Pa_s_m3"]
    re_per_q = p["rho"]*p["pipe_diameter"]/(p["water_viscosity"]*h["pipe_area_m2"])
    re = re_per_q*qabs
    slope = p["Rh"]
    if h["nonlinear_pipe"] and re > h["Re_laminar"]:
        span = h["Re_turbulent"]-h["Re_laminar"]
        s = min((re-h["Re_laminar"])/span, 1.)
        w = s*s*(3-2*s); dw_dre = 6*s*(1-s)/span
        arg = 6.9/re + (h["pipe_absolute_roughness_m"]/(3.7*p["pipe_diameter"]))**1.11
        denominator = -1.8*np.log10(arg)
        f = denominator**-2
        df_dre = -2*denominator**-3 * (1.8*6.9/(np.log(10)*arg*re**2))
        C = h["pipe_length_m"]/p["pipe_diameter"]*p["rho"]/(2*h["pipe_area_m2"]**2)
        turb = C*f*qabs*qabs
        dturb = C*(2*f*qabs+df_dre*re_per_q*qabs*qabs)
        slope += dw_dre*re_per_q*(turb-Rlam*qabs)+w*(dturb-Rlam)
    return slope+2*h["additional_component_K_Pa_s2_m6"]*qabs


def outputs(t, z, p, h, *, force_N=None, actuator=None, segment=None):
    x, v, V, Q = np.asarray(z)[:4]
    if actuator is not None:
        if force_N is not None or segment is None:
            raise ValueError("Use exactly one input mode: force fixture or STEP/DIR.")
        o = step4.outputs(t, np.asarray(z)[:3], p, actuator, segment)
        M = step4.total_equiv_mass(actuator)
    else:
        if force_N is None:
            raise ValueError("A force input is required in fixture mode.")
        # A prescribed force profile is independent of the state. The inherited
        # one-second loading ramp still acts only once, at the global start.
        force_value = force_N(t) if callable(force_N) else force_N
        if not np.all(np.isfinite(force_value)) or np.any(np.asarray(force_value) < 0):
            raise ValueError("Force profile must be finite and nonnegative.")
        o = step4.legacy_plant.outputs(t, np.asarray(z)[:3], force_value, p)
        o["pressure"] = o["p"]
        M = p["mass_equivalent"]  # Do not insert the motor into a hanging-load fixture.
    loss = hydraulic_loss(Q, p, h)
    acceleration = (o["Fm"]-p["damping"]*v-o["Ffriction"]-o["Fn"])/M
    drive = o["pressure"] + p["rho"]*p["gravity"]*p["head"]
    Qdot = (drive-loss["dp_loss"])/h["Ih_Pa_s2_m3"]
    o.update(loss)
    o.update(x=x, v=v, V=V, Q=Q, Qdot=Qdot, acceleration=acceleration,
             drive=drive, Q_old_algebraic=np.maximum(drive, 0.)/p["Rh"],
             ydot=Q/o["A"], fluid_kinetic_energy=.5*h["Ih_Pa_s2_m3"]*Q**2,
             hydraulic_dissipation_W=loss["dp_loss"]*Q,
             Dcage=p["D0"]-x/np.pi, Lring=p["L0"]-x,
             radial_gap=(o["y"]-x)/(2*np.pi))
    return o


def rhs(t, z, p, h, *, force_N=None, actuator=None, segment=None):
    o = outputs(t, z, p, h, force_N=force_N, actuator=actuator, segment=segment)
    dz = [o["v"], o["acceleration"], -o["Q"], o["Qdot"]]
    if len(z) == 6:
        power = (o["Fm"]*o["v"] + p["rho"]*p["gravity"]*p["head"]*o["Q"]
                 - p["damping"]*o["v"]**2 - o["Ffriction"]*o["v"]
                 - o["hydraulic_dissipation_W"])
        dz += [o["Q"], power]
    return np.asarray(dz)


def jacobian(t, z, p, h, *, force_N=None, actuator=None, segment=None):
    o = outputs(t, z, p, h, force_N=force_N, actuator=actuator, segment=segment)
    x, v, V, Q = z[:4]
    area = float(o["A"])
    contact_k = p["contact_stiffness"] if o["overlap"] > 0 else 0.
    fnV = contact_k/area
    darea = 2/(np.pi*o["D"])
    dpdx = contact_k/area
    dpdV = (step4.legacy_plant.pressure_derivative(V,p)
                                  + fnV/area-o["Fn"]*darea/area**2)
    dFm_dx = dFm_dv = 0.
    M = p["mass_equivalent"] if actuator is None else step4.total_equiv_mass(actuator)
    if actuator is not None and segment.enabled:
        radius = actuator["radius_m"]
        dFm_dx = -o["torque_cap"]*np.cos(o["phase_error"])*(actuator["full_steps_per_rev"]/4)/radius**2
        w = abs(v/radius)
        speeds, caps = np.asarray(actuator["torque_speed_rad_s"]), np.asarray(actuator["torque_Nm"])
        cap_slope = 0.
        if speeds[0] <= w < speeds[-1]:
            j = min(np.searchsorted(speeds,w,side="right")-1,len(speeds)-2)
            cap_slope = (caps[j+1]-caps[j])/(speeds[j+1]-speeds[j])
        dFm_dv = cap_slope*np.sign(v)*np.sin(o["phase_error"])/radius**2
    dfr = p["friction"]/p["friction_velocity"]*(1-np.tanh(v/p["friction_velocity"])**2)
    dloss = hydraulic_slope(Q,p,h)
    J = np.zeros((len(z),len(z)))
    J[0,1] = 1.
    J[1,:3] = [(dFm_dx-contact_k)/M,(dFm_dv-p["damping"]-dfr)/M,-fnV/M]
    J[2,3] = -1.
    J[3,[0,2,3]] = [dpdx/h["Ih_Pa_s2_m3"],dpdV/h["Ih_Pa_s2_m3"],-dloss/h["Ih_Pa_s2_m3"]]
    if len(z)==6:
        J[4,3] = 1.
        J[5,[0,1,3]] = [v*dFm_dx, o["Fm"]+v*dFm_dv-2*p["damping"]*v-o["Ffriction"]-v*dfr,
                         p["rho"]*p["gravity"]*p["head"]-o["dp_loss"]-Q*dloss]
    return J


def _guards(t, z, p, h, *, force_N=None, actuator=None, segment=None, target_mL=None):
    if actuator is None:
        names = ["50mm_cage_limit", "sphere_volume_limit", "extension_outside_model"]
        values = [p["xmax"]-z[0], z[2]-p["Vmin"], z[0]+1e-6]
    else:
        names = list(step4.STOP_NAMES)
        values = list(step4.guard_values(t, np.asarray(z)[:3], p, actuator, segment))
    names.append("forward_liquid_flow_boundary")
    values.append(z[3])
    if target_mL is not None:
        names.append("target_discharge_reached")
        values.append(target_mL*1e-6 - (p["V0"]-z[2]))
    return names, np.asarray(values)


def _simulate(p, h, schedules, *, force_N=None, actuator=None, initial=None,
              max_step=.25, target_mL=None):
    z = np.array([0., 0., p["V0"], 0., 0., 0.] if initial is None else initial, dtype=float)
    if z.shape != (6,) or not np.all(np.isfinite(z)) or z[3] < 0:
        raise ValueError("Supply x,v,V,Q plus two audit integrals; initial Q must be nonnegative.")
    if target_mL is not None and not 0 < target_mL < 1e6*(p["V0"]-p["Vmin"]):
        raise ValueError("Target must be positive and below the geometric sphere limit.")
    pieces, traces = [], []
    stop = "time_horizon"
    for start, end, segment in schedules:
        if segment is not None and not segment.valve_open:
            stop = "valve_closure_not_modelled"
            break  # Do not reset Q or silently destroy fluid kinetic energy.
        kw = dict(force_N=force_N, actuator=actuator, segment=segment)
        names, values = _guards(start, z, p, h, **kw, target_mL=target_mL)
        for name, value in zip(names, values):
            if name != "forward_liquid_flow_boundary" and value <= 0:
                stop = name
                break
        if stop != "time_horizon":
            break
        if z[3] <= 0 and outputs(start, z, p, h, **kw)["drive"] <= 0:
            stop = "initial_forward_flow_not_available"
            break
        events = []
        for j in range(len(names)):
            def event(t, state, j=j):
                return _guards(t, state, p, h, **kw, target_mL=target_mL)[1][j]
            event.terminal = True
            event.direction = -1
            events.append(event)
        sol = solve_ivp(lambda t, state: rhs(t, state, p, h, **kw), (start, end), z,
                        method="Radau", rtol=2e-8,
                        atol=[2e-11, 2e-10, 2e-13, 2e-14, 2e-13, 2e-10],
                        max_step=max_step, dense_output=True, events=events,
                        jac=lambda t,state: jacobian(t,state,p,h,**kw))
        if not sol.success:
            raise RuntimeError(sol.message)
        t = np.unique(np.r_[sol.t, np.linspace(sol.t[0], sol.t[-1],
                                              max(2, int((sol.t[-1]-sol.t[0])/.1)+1))])
        states = sol.sol(t)
        out = outputs(t, states, p, h, **kw)
        out.update(t=t, z=states)
        traces.append(out)
        pieces.append(sol)
        z = states[:, -1]
        for name, times in zip(names, sol.t_events):
            if len(times):
                stop = name
                break
        if stop != "time_horizon":
            break
    if not traces:
        raise ValueError("Initial state/schedule is outside the reference scope: " + stop)
    data = {key: np.concatenate([tr[key] for tr in traces], axis=1 if key == "z" else 0)
            for key in traces[0]}
    keep = np.r_[np.diff(data["t"]) > 1e-12, True]
    for key in data:
        data[key] = data[key][:, keep] if key == "z" else data[key][keep]
    M = p["mass_equivalent"] if actuator is None else step4.total_equiv_mass(actuator)
    E = (.5*M*data["v"]**2 + .5*p["contact_stiffness"]*data["overlap"]**2
         + step4.legacy_baseline.shell_energy(data["V"], p) + data["Ug"]
         + data["fluid_kinetic_energy"])
    data["energy_J"] = E
    report = dict(stop=stop, end_s=float(data["t"][-1]),
                  input_mode="force_fixture" if actuator is None else "step_dir",
                  physical_states=["x_m", "v_m_s", "V_m3", "Q_m3_s"],
                  numerical_audit_states=2, physical_hardware_validation=False,
                  illustrative_actuator=bool(actuator and actuator.get("illustrative_only", False)),
                  material_relaxation_identified=False, component_nonlinearity_identified=False,
                  original_Rh_preserved=float(p["Rh"]), total_equiv_mass_kg=float(M),
                  min_flow_mL_s=float(np.min(data["Q"])*1e6), max_Re=float(np.max(data["Re"])),
                  max_hydraulic_dissipation_W=float(np.max(data["hydraulic_dissipation_W"])),
                  discharged_mL=float((data["V"][0]-data["V"][-1])*1e6),
                  x_mm=float(data["x"][-1]*1000),
                  below_natural_volume=bool(np.any(data["V"] < p["Vnatural"])),
                  mass_balance_error_mL=float(np.max(np.abs(data["z"][4]-data["z"][4,0]
                                               -(data["V"][0]-data["V"])))*1e6),
                  energy_balance_error_J=float(np.max(np.abs(E-E[0]-(data["z"][5]-data["z"][5,0])))),
                  finite=bool(all(np.all(np.isfinite(v)) for v in data.values())),
                  post_boundary_motion_predicted=False, valve_closure_transient_predicted=False)
    if actuator:
        report.update(max_abs_phase_error_rad=float(np.max(np.abs(data["phase_error"]))),
                      min_required_cable_tension_N=float(np.min(data["cable_tension"])))
    return dict(data=data, report=report, solutions=pieces,
                segments=[asdict(s) for _, _, s in schedules if s is not None])


def simulate_force(force_N, horizon_s, p, h, *, initial=None, max_step=2., target_mL=None):
    if horizon_s <= 0 or not np.isfinite(force_N) or force_N < 0:
        raise ValueError("Use a nonnegative finite force and positive horizon.")
    return _simulate(p, h, [(0., horizon_s, None)], force_N=force_N, initial=initial,
                     max_step=max_step, target_mL=target_mL)


def simulate_stepper(commands, horizon_s, p, h, actuator, *, initial=None,
                     max_step=.025, target_mL=None):
    step4.validate_actuator(actuator)
    segments = step4.shape_commands(commands, horizon_s, actuator)
    return _simulate(p, h, [(s.start_s, s.end_s, s) for s in segments], actuator=actuator,
                     initial=initial, max_step=max_step, target_mL=target_mL)


def simulate_force_profile(profile, horizon_s, p, h, *, force_limit_N,
                           initial=None, max_step=.25, target_mL=None):
    """Synthetic force-port experiment, not a measured motor torque envelope.

    The supplied profile must support scalar and vector time. The force cap is
    checked at every RHS/output evaluation, not only at output sample times.
    """
    if not callable(profile) or not np.isfinite(horizon_s) or horizon_s <= 0:
        raise ValueError("A callable force profile and finite positive horizon are required.")
    if not np.isfinite(force_limit_N) or force_limit_N <= 0:
        raise ValueError("The scenario force limit must be finite and positive.")

    def checked(t):
        value = np.asarray(profile(t), dtype=float)
        if value.shape != np.asarray(t).shape:
            raise ValueError("Force profile must preserve scalar/vector time shape.")
        if not np.all(np.isfinite(value)) or np.any(value < 0) or np.any(value > force_limit_N + 1e-12):
            raise ValueError("Force profile exceeds its declared nonnegative bound.")
        return value

    checked(0.0)
    checked(np.array([0.0, horizon_s]))
    result = _simulate(p, h, [(0., horizon_s, None)], force_N=checked,
                       initial=initial, max_step=max_step, target_mL=target_mL)
    result["report"].update(input_mode="prescribed_bounded_force_profile",
                            force_limit_N=float(force_limit_N),
                            hardware_force_limit_identified=False,
                            data_kind="synthetic_pilot_not_experimental")
    return result


def sample_state(result, times):
    t = np.atleast_1d(times).astype(float)
    if t.size == 0 or not np.all(np.isfinite(t)) or np.min(t) < 0 or np.max(t) > result["report"]["end_s"]+1e-8:
        raise ValueError("Do not extrapolate beyond a model-validity stop.")
    z = np.empty((6, len(t)))
    for sol in result["solutions"]:
        use = (t >= sol.t[0]-1e-10) & (t <= sol.t[-1]+1e-10)
        if np.any(use):
            z[:, use] = sol.sol(t[use])
    return z
