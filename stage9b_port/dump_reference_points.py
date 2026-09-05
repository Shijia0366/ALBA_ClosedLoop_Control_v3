"""Dump reference right-hand-side / Jacobian / observable values from the FROZEN
Stage-8 five-state reference model, for pointwise verification of the MATLAB port.

Nothing here modifies the frozen model. It constructs stop_study.Plant exactly as
Stage 8 does (which is what applies the +0.10 m head and the conditional torque
envelope over the illustrative actuator file) and evaluates it on a state grid.

The MATLAB port must reproduce these numbers; it is not allowed to be "close in
spirit". The tolerance used downstream is the same 1e-8 relative bound the
Stage-8 package uses for its own regression_checks.

Usage:  python dump_reference_points.py <path to ALBA_Step8_SINDy> <out.json>
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np


class FixedCommand:
    """Minimal profile: the plant only asks for (count, rate) at a time.

    Using a constant command makes the comparison a pure function of the state,
    so a mismatch cannot be blamed on a differently interpolated input spline.
    """

    enabled = True
    valve_open = True
    requested_steps_s = 0.0

    def __init__(self, count, rate):
        self._count = float(count)
        self._rate = float(rate)

    def scalar(self, t):
        return self._count, self._rate

    def pulse_state(self, t):
        return self._count, self._rate


def main():
    pkg = Path(sys.argv[1])
    out = Path(sys.argv[2])
    sys.path.insert(0, str(pkg / "frozen_reference"))
    import stop_study as ss  # noqa: E402

    plant = ss.Plant()
    p, h, a, cfg = plant.p, plant.h, plant.a, plant.cfg
    r = cfg.radius_m
    spr = cfg.steps_per_rev

    meta = dict(
        total_equiv_mass_kg=float(plant.M),
        rotor_inertia_kg_m2=float(a["rotor_inertia_kg_m2"]),
        nonmotor_equiv_mass_kg=float(a["nonmotor_equiv_mass_kg"]),
        torque_cap_Nm=float(cfg.torque_cap),
        radius_m=float(r),
        steps_per_rev=int(spr),
        full_steps_per_rev=int(a["full_steps_per_rev"]),
        microsteps_per_full_step=int(a["microsteps_per_full_step"]),
        supply_voltage_V=float(a["supply_voltage_V"]),
        current_assumed_A=float(a["current_assumed_A"]),
        current_limit_measured=bool(a["current_limit_measured"]),
        torque_curve_matches_drive_conditions=bool(a["torque_curve_matches_drive_conditions"]),
        damping=float(p["damping"]),
        friction=float(p["friction"]),
        friction_velocity=float(p["friction_velocity"]),
        contact_stiffness=float(p["contact_stiffness"]),
        head_Pa=float(plant.head_Pa),
        Ih=float(plant.I),
        Rh=float(p["Rh"]),
        xmax=float(p["xmax"]),
        Vmin=float(p["Vmin"]),
        D0=float(p["D0"]),
        Vsilicone=float(p["Vsilicone"]),
        legacy_force_max_N_DO_NOT_USE=float(p["force_max"]),
        max_rate=float(cfg.max_rate),
        max_accel=float(cfg.max_accel),
        closed_brake_accel=float(cfg.closed_brake_accel),
        note=(
            "torque_cap overrides the illustrative 0.0375 N*m file value. It is a "
            "conditional flat envelope, NOT a measured 12 V running-torque curve. "
            "force_max is a legacy fixture limit and must not be used as a motor bound."
        ),
    )

    rows = []
    volumes_mL = [500.0, 420.0, 340.0, 260.0, 180.0, 110.0]
    flows_mLs = [0.0, 2.0, 8.0, 15.0, 20.0]
    speeds = [-0.004, 0.0, 0.0015, 0.012]
    overlaps = [0.0, 1e-5, 1.2e-4, 5.5e-4]
    phases = [-1.2, -0.4, 0.0, 0.45, 1.25]

    for VmL in volumes_mL:
        V = VmL * 1e-6
        D, A, pe, pg, y, dp = plant.geometry(V)
        for QmLs in flows_mLs:
            Q = QmLs * 1e-6
            for v in speeds:
                for ov in overlaps:
                    x = y + ov
                    for phi in phases:
                        count = (phi / 50.0 + x / r) * spr / (2 * math.pi)
                        prof = FixedCommand(count, 0.0)
                        z = np.array([x, v, V, Q, 0.0, 0.0])
                        for sealed in (False, True):
                            o = plant.obs(8.0, z, prof, sealed=sealed)
                            f = plant.rhs(8.0, z, prof, sealed=sealed)
                            J = plant.jac(8.0, z, prof, sealed=sealed)
                            g = plant.guard_values(8.0, z, prof, sealed=sealed)
                            rows.append(dict(
                                x=x, v=v, V=V, Q=Q, count=count, sealed=int(sealed),
                                phi=float(o["phi"]), torque=float(o["torque"]),
                                Fm=float(o["Fm"]), Fn=float(o["Fn"]),
                                fr=float(o["fr"]), acc=float(o["acc"]),
                                tension=float(o["tension"]), loss=float(o["loss"]),
                                slope=float(o["slope"]), Qdot=float(o["Qdot"]),
                                pe=float(o["pe"]), pg=float(o["pg"]),
                                A=float(o["A"]), D=float(o["D"]), y=float(o["y"]),
                                rhs=[float(q) for q in f],
                                jac=[[float(q) for q in rowJ] for rowJ in J],
                                guards=[float(q) for q in g],
                            ))

    # a static loss sweep as well, since the waterway is shared with the ROM
    loss_rows = []
    for QmLs in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0]:
        Q = QmLs * 1e-6
        L, sl = plant.loss(Q)
        loss_rows.append(dict(Q_m3_s=Q, loss_Pa=float(L), slope=float(sl)))

    out.write_text(json.dumps(dict(meta=meta, points=rows, loss=loss_rows), indent=1))
    print(f"wrote {out} : {len(rows)} state points, {len(loss_rows)} loss points")


if __name__ == "__main__":
    main()
