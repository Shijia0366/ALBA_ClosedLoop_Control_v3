"""ALBA STEP/DIR actuator, coupled to the unchanged spherical bladder plant.

This is a phase-averaged, first-harmonic, taut-transmission approximation.
It is NOT a PWM/current-regulator simulation or a validated lost-step model.
Stop at the +/-pi/2 electrical-phase boundary; do not fabricate post-slip motion.
The physical states remain x, v, V. Command phase is a known integral of the
shaped pulse rate; two optional integrals audit discharged volume and work.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
import json
import sys

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "plant_snapshot"))
import baseline as legacy_baseline
import selfweight_model as legacy_plant


class MissingHardwareParameters(ValueError):
    pass


def load_plant():
    p = json.loads((ROOT / "plant_parameters.json").read_text(encoding="utf-8"))
    if not np.isclose(p["Rh"], 1644008304.5256217, rtol=1e-12):
        raise ValueError("Unexpected hydraulic baseline; do not silently refit Rh.")
    return p


def load_actuator(name="hardware_pending.json", *, allow_illustrative=False):
    a = json.loads((ROOT / name).read_text(encoding="utf-8"))
    if a.get("illustrative_only", False):
        if not allow_illustrative:
            raise MissingHardwareParameters("Illustrative parameters require explicit opt-in.")
    else:
        missing = [k for k in ("microsteps_per_full_step", "driver_board",
                   "supply_voltage_V", "current_limit_A", "torque_speed_rad_s", "torque_Nm")
                   if a.get(k) is None]
        missing += [k for k in ("current_motor_confirmed", "flashed_firmware_confirmed",
                               "torque_curve_matches_drive_conditions") if a.get(k) is not True]
        if missing:
            raise MissingHardwareParameters("Please confirm: " + ", ".join(missing))
    validate_actuator(a)
    return a


def validate_actuator(a):
    for k in ("full_steps_per_rev", "microsteps_per_full_step", "radius_m",
              "rotor_inertia_kg_m2", "max_step_rate_s", "acceleration_steps_s2"):
        if not np.isfinite(a[k]) or a[k] <= 0:
            raise ValueError(f"Invalid actuator parameter {k}")
    if a["full_steps_per_rev"] % 4 != 0:
        raise ValueError("This two-phase approximation requires four full steps per electrical cycle.")
    if a["microsteps_per_full_step"] not in (1, 2, 4, 8, 16):
        raise ValueError("Unsupported A4988 subdivision.")
    if a["nonmotor_equiv_mass_kg"] < 0:
        raise ValueError("Nonmotor reflected inertia must not be negative.")
    speeds, torques = np.asarray(a["torque_speed_rad_s"]), np.asarray(a["torque_Nm"])
    if (speeds.size < 2 or speeds.shape != torques.shape or speeds[0] != 0 or
            np.any(np.diff(speeds) <= 0) or np.any(torques <= 0) or
            not np.all(np.isfinite(speeds)) or not np.all(np.isfinite(torques))):
        raise ValueError("Supply a finite positive torque envelope starting at zero speed.")


@dataclass(frozen=True)
class Command:
    time_s: float
    requested_steps_s: float
    enabled: bool = True
    valve_open: bool = True


@dataclass(frozen=True)
class Segment:
    start_s: float
    end_s: float
    rate0_steps_s: float
    slope_steps_s2: float
    pulses0: float
    requested_steps_s: float
    limited_steps_s: float
    enabled: bool
    valve_open: bool

    def pulse_state(self, t):
        dt = np.clip(np.asarray(t) - self.start_s, 0., self.end_s - self.start_s)
        rate = self.rate0_steps_s + self.slope_steps_s2 * dt
        count = self.pulses0 + self.rate0_steps_s * dt + .5 * self.slope_steps_s2 * dt**2
        return count, rate


def shape_commands(commands, horizon_s, actuator):
    """Continuous-time analogue of the inspected firmware ramp.

    Zero/disable stops pulse generation immediately; a sign change ramps via
    zero. Timer quantization, the ~50 us DIR dwell and service faults are not
    reconstructed. An applied rate is NOT an encoder-measured speed.
    """
    if not commands or commands[0].time_s != 0 or horizon_s <= commands[-1].time_s:
        raise ValueError("Commands must start at zero and precede the horizon.")
    times = [c.time_s for c in commands]
    if np.any(np.diff(times) <= 0):
        raise ValueError("Command times must increase strictly.")
    segments = []
    rate = count = 0.
    accel = actuator["acceleration_steps_s2"]
    cap = actuator["max_step_rate_s"]
    for i, cmd in enumerate(commands):
        if not np.isfinite(cmd.requested_steps_s):
            raise ValueError("Finite pulse commands are required.")
        t = float(cmd.time_s)
        end = float(commands[i + 1].time_s if i + 1 < len(commands) else horizon_s)
        target = float(np.clip(cmd.requested_steps_s, -cap, cap))
        if not cmd.enabled or target == 0.:
            rate = 0.
            segments.append(Segment(t, end, 0., 0., count, cmd.requested_steps_s,
                                    target, cmd.enabled, cmd.valve_open))
            continue
        while t < end - 1e-12:
            goal = 0. if rate * target < 0 else target
            if abs(rate - goal) < 1e-10:
                rate = goal
                slope, dt = 0., end - t
            else:
                slope = float(np.sign(goal - rate) * accel)
                dt = min(end - t, abs(goal - rate) / accel)
            segments.append(Segment(t, t + dt, rate, slope, count,
                                    cmd.requested_steps_s, target, cmd.enabled, cmd.valve_open))
            count += rate * dt + .5 * slope * dt**2
            rate += slope * dt
            t += dt
            if abs(rate) < 1e-10:
                rate = 0.
    return segments


def total_equiv_mass(a):
    # Nonmotor load is specified separately; the old lumped mass is not
    # silently treated as a measured total and then counted a second time.
    return a["nonmotor_equiv_mass_kg"] + a["rotor_inertia_kg_m2"] / a["radius_m"]**2


def outputs(t, z, p, a, segment, *, connected=True):
    x, v, V = np.asarray(z)[:3]
    g = legacy_plant.geometry(V, p)
    count, pulse_rate = segment.pulse_state(t)
    pulses_per_rev = a["full_steps_per_rev"] * a["microsteps_per_full_step"]
    theta_cmd = 2 * np.pi * count / pulses_per_rev
    r = a["radius_m"]
    omega = v / r
    phi = (a["full_steps_per_rev"] / 4) * (theta_cmd - x / r)
    torque_cap = np.interp(np.abs(omega), a["torque_speed_rad_s"], a["torque_Nm"])
    torque = torque_cap * np.sin(phi) if segment.enabled else np.zeros_like(phi)
    motor_force = torque / r
    overlap = np.maximum(x - g["y"], 0.) if connected else np.zeros_like(x)
    fn = p["contact_stiffness"] * overlap
    pressure = g["pe"] + g["pg"] + fn / g["A"]
    drive = pressure + p["rho"] * p["gravity"] * p["head"]
    Q = np.maximum(drive, 0.) / p["Rh"] if segment.valve_open and connected else np.zeros_like(V)
    fr = p["friction"] * np.tanh(v / p["friction_velocity"])
    acceleration = (motor_force - p["damping"] * v - fr - fn) / total_equiv_mass(a)
    # In this reduced transmission, baseline friction/damping are assigned
    # to the nonmotor fixture. A cable cannot provide negative tension.
    cable_tension = a["nonmotor_equiv_mass_kg"] * acceleration + p["damping"] * v + fr + fn
    g.update(x=x, v=v, V=V, Fn=fn, overlap=overlap, pressure=pressure,
             Q=Q, Fm=motor_force, Ffriction=fr, torque=torque, torque_cap=torque_cap,
             omega=omega, phase_error=phi, command_angle=theta_cmd, command_pulses=count,
             requested_steps_s=np.zeros_like(x) + segment.requested_steps_s,
             applied_steps_s=np.zeros_like(x) + pulse_rate,
             commanded_cable_speed=r * 2 * np.pi * pulse_rate / pulses_per_rev,
             motor_pos_rev=x / (2 * np.pi * r), motor_vel_rev_s=omega / (2 * np.pi),
             enabled=np.zeros_like(x) + segment.enabled,
             valve_open=np.zeros_like(x) + segment.valve_open,
             cable_tension=cable_tension)
    return g


def rhs(t, z, p, a, segment, *, connected=True):
    o = outputs(t, z, p, a, segment, connected=connected)
    dz = [o["v"], (o["Fm"] - p["damping"] * o["v"] - o["Ffriction"] - o["Fn"])
          / total_equiv_mass(a), -o["Q"]]
    if len(z) == 5:
        power = (o["Fm"] * o["v"] + p["rho"] * p["gravity"] * p["head"] * o["Q"]
                 - p["damping"] * o["v"]**2 - o["Ffriction"] * o["v"] - p["Rh"] * o["Q"]**2)
        dz += [o["Q"], power]
    return np.asarray(dz)


STOP_NAMES = ("synchronism_margin_exhausted", "50mm_cage_limit", "sphere_volume_limit",
              "extension_outside_model", "torque_envelope_speed_range_exceeded",
              "taut_cable_model_boundary")


def guard_values(t, z, p, a, segment, connected=True):
    o = outputs(t, z, p, a, segment, connected=connected)
    return np.asarray([np.pi / 2 - abs(o["phase_error"]) if segment.enabled else np.pi / 2,
                       p["xmax"] - z[0], z[2] - p["Vmin"], z[0] + 1e-6,
                       (a["torque_speed_rad_s"][-1] - abs(z[1] / a["radius_m"]))
                       if segment.enabled else a["torque_speed_rad_s"][-1],
                       o["cable_tension"] + 1e-5])  # numerical tolerance, not a physical preload


def simulate(commands, horizon_s, p, a, *, connected=True, initial=None, max_step=.025):
    validate_actuator(a)
    segments = shape_commands(commands, horizon_s, a)
    z = np.array([0., 0., p["V0"], 0., 0.] if initial is None else initial, dtype=float)
    if z.shape != (5,):
        raise ValueError("Initial x,v,V and two zero audit integrals are required.")
    traces = []
    stop = "time_horizon"
    for segment in segments:
        guards = guard_values(segment.start_s, z, p, a, segment, connected)
        bad = np.flatnonzero(guards <= 0)
        if bad.size:
            stop = STOP_NAMES[bad[0]]
            break
        events = []
        for j in range(len(STOP_NAMES)):
            def event(t, state, j=j):
                return guard_values(t, state, p, a, segment, connected)[j]
            event.terminal = True
            event.direction = -1
            events.append(event)
        sol = solve_ivp(lambda t, state: rhs(t, state, p, a, segment, connected=connected),
                        (segment.start_s, segment.end_s), z, method="Radau",
                        rtol=2e-8, atol=[2e-11, 2e-10, 2e-13, 2e-13, 2e-10],
                        max_step=max_step, dense_output=True, events=events)
        if not sol.success:
            raise RuntimeError(sol.message)
        t = np.unique(np.r_[sol.t, np.linspace(sol.t[0], sol.t[-1], max(2, int((sol.t[-1]-sol.t[0])/.02)+1))])
        states = sol.sol(t)
        o = outputs(t, states, p, a, segment, connected=connected)
        o.update(t=t, z=states)
        traces.append(o)
        z = states[:, -1]
        for name, times in zip(STOP_NAMES, sol.t_events):
            if len(times):
                stop = name
                break
        if stop != "time_horizon":
            break
    if not traces:
        raise ValueError("Initial condition or command transition is outside model validity: " + stop)
    data = {key: np.concatenate([tr[key] for tr in traces], axis=1 if key == "z" else 0)
            for key in traces[0]}
    # Remove boundary duplicates, retaining the new command/valve state.
    keep = np.r_[np.diff(data["t"]) > 1e-12, True]
    for key in data:
        data[key] = data[key][:, keep] if key == "z" else data[key][keep]
    E = .5 * total_equiv_mass(a) * data["v"]**2 + .5 * p["contact_stiffness"] * data["overlap"]**2
    if connected:
        E += legacy_baseline.shell_energy(data["V"], p) + data["Ug"]
    mass_error = np.max(np.abs(data["z"][3] - (data["V"][0] - data["V"]))) * 1e6
    energy_error = np.max(np.abs(E - E[0] - (data["z"][4] - data["z"][4, 0])))
    report = {
        "stop": stop, "end_s": float(data["t"][-1]),
        "illustrative_only": bool(a.get("illustrative_only", False)),
        "physical_hardware_validation": False,
        "post_synchronism_loss_predicted": False,
        "total_equiv_mass_kg": total_equiv_mass(a),
        "mass_balance_error_mL": float(mass_error),
        "mechanical_hydraulic_energy_error_J": float(energy_error),
        "max_abs_phase_error_rad": float(np.max(np.abs(data["phase_error"]))),
        "max_abs_torque_Nm": float(np.max(np.abs(data["torque"]))),
        "max_torque_bound_violation_Nm": float(np.max(np.maximum(np.abs(data["torque"]) - data["torque_cap"], 0.))),
        "min_required_cable_tension_N": float(np.min(data["cable_tension"])),
        "post_slack_motion_predicted": False,
        "discharged_mL": float((data["V"][0] - data["V"][-1]) * 1e6),
        "x_mm": float(data["x"][-1] * 1000),
        "endpoint_inside_fixed_id_25_to_250mL_window": bool(25e-6 <= (p["V0"] - data["V"][-1]) <= 250e-6),
        "finite": bool(all(np.all(np.isfinite(v)) for v in data.values())),
    }
    return {"data": data, "report": report, "segments": [asdict(s) for s in segments]}
