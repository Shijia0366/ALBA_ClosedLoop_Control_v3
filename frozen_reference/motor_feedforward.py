"""Inverse-dynamics motor command; the forward plant still solves actual flow.

This is a nominal model-based feedforward design, not experimentally fitted
control or an ideal imposed flow source. The contact load starts at zero and
becomes positive immediately. No measured or invented initial preload is used.
All STEP commands are continuous, averaged microstep-phase commands.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline


def smooth_weight(t, start, duration):
    u = min(max((t-start)/duration, 0.0), 1.0)
    if u <= 0.0:
        return 0.0, 0.0, 0.0
    if u >= 1.0:
        return 1.0, 0.0, 0.0
    return (u**3*(10-15*u+6*u*u),
            30*u*u*(1-u)**2/duration,
            60*u*(1-u)*(1-2*u)/duration**2)


def startup_weight(t, duration):
    """Zero initial load/rate, finite initial acceleration, C2 end join."""
    u = min(max(t/duration, 0.0), 1.0)
    if u >= 1.0:
        return 1.0, 0.0, 0.0
    return (6*u*u-8*u**3+3*u**4,
            12*u*(1-u)**2/duration,
            12*(1-u)*(1-3*u)/duration**2)


class LoadAwareProfile:
    """Plan compressive contact loading and convert it to a STEP phase input.

    The two-state planning ODE is used only to generate a command. A separate
    six-state plant must replay that command, check contact, and audit energy.
    The natural open-valve discharge is not mistaken for lost cable motion.
    """
    enabled = True
    valve_open = True
    requested_steps_s = 0.0

    def __init__(self, plant, *, base=None, brake_start=None):
        self.plant, self.cfg, self.p = plant, plant.cfg, plant.p
        self.base, self.brake_start = base, brake_start
        c, p = self.cfg, self.p
        self.m_per_step = 2*math.pi*c.radius_m/c.steps_per_rev
        self.full_loss = plant.loss(c.q_target)[0]
        self.tail_loss = plant.loss(c.q_close_design)[0]
        if brake_start is None:
            start, stop = 0.0, 60.0
            initial = [p['V0'], 0.0]
        else:
            if base is None or brake_start < c.rise_s:
                raise ValueError('Slowdown requires the completed startup profile')
            start = float(brake_start)
            stop = start+c.slowdown_s+4.0
            initial = base.plan.sol(start)

        def rhs(t, state):
            V, Q = state
            geom = plant.geometry(V)
            F, _, _ = self.load(t, V, Q, None, geom)
            H = geom[2]+geom[3]+plant.head_Pa
            return [-Q, (H+F/geom[1]-plant.loss(Q)[0])/plant.I]

        def domain(t, state):
            return state[0]-p['Vmin']-1e-9
        domain.terminal = True
        domain.direction = -1

        def jac(t, state):
            # During the fully active loading/slowdown plan, H(V) cancels
            # analytically. A finite-difference Jacobian can chase roundoff
            # in that zero derivative with unphysically large V perturbations.
            s = 1.0 if brake_start is not None else startup_weight(t, c.rise_s)[0]
            Hp = plant.geometry(state[0])[5]
            return [[0., -1.], [(1-s)*Hp/plant.I, -plant.loss(state[1])[1]/plant.I]]

        self.plan = solve_ivp(rhs, (start, stop), initial, method='Radau',
                              rtol=3e-10, atol=[2e-14, 1e-15], max_step=.02,
                              dense_output=True, events=domain, jac=jac)
        if not self.plan.success:
            raise RuntimeError(self.plan.message)
        end = self.plan.t[-1]
        # Resolve the hydraulic inertial transient and the regularised
        # friction near zero velocity, not just the seconds-long loading ramp.
        early = np.r_[np.arange(0., .02, .0001),
                      np.arange(.02, .1, .0002), np.arange(.1, 1., .002)]
        grid = np.unique(np.r_[np.arange(start, end, .005),
                               early[(early >= start)&(early <= end)], start, end])
        grid = grid[np.r_[np.diff(grid) > 1e-9, True]]
        samples = [self.reference(float(t)) for t in grid]
        raw_count = np.array([r['count'] for r in samples])
        count = np.maximum.accumulate(raw_count)
        self.monotone_adjustment_steps = float(max(count-raw_count))
        boundary = ((1, 0.0), 'not-a-knot') if start == 0 else 'not-a-knot'
        self.command = CubicSpline(grid, count, bc_type=boundary, extrapolate=False)
        self.grid = grid
        self.reference_peak_torque = max(abs(r['torque']) for r in samples)

    def load(self, t, V, Q, Qdot, geom=None):
        """Return F, dF/dt, d2F/dt2 along the planning trajectory."""
        p, plant, c = self.p, self.plant, self.cfg
        D, A, pe, pg, _, Hp = plant.geometry(V) if geom is None else geom
        H = pe+pg+plant.head_Pa
        W = V+p['Vsilicone']
        ua = p['Dref_inner']/(2*(3*V/(4*math.pi))**(1/3))
        ub = p['Dref_outer']/D
        Hpp = 2*p['mu']/9*(-(4*ua+7*ua**4)/V**2+(4*ub+7*ub**4)/W**2)
        Hpp += 10*pg/(9*W**2)
        Ap, App = 2/(math.pi*D), -4/(math.pi**2*D**4)
        if self.brake_start is not None and t >= self.brake_start:
            s, ds, dds = smooth_weight(t, self.brake_start, c.slowdown_s)
            weights = [(1-s, -ds, -dds, self.full_loss),
                       (s, ds, dds, self.tail_loss)]
        else:
            s, ds, dds = startup_weight(t, c.rise_s)
            weights = [(s, ds, dds, self.full_loss)]
        F = Fdot = Fddot = 0.0
        for w, dw, ddw, L in weights:
            C = A*(L-H)
            Cp = Ap*(L-H)-A*Hp
            Cpp = App*(L-H)-2*Ap*Hp-A*Hpp
            F += w*C
            Fdot += dw*C-w*Cp*Q
            if Qdot is not None:
                Fddot += ddw*C-2*dw*Cp*Q+w*(Cpp*Q*Q-Cp*Qdot)
        if F < -1e-8:
            raise ValueError('Requested pressure ramp would require tensile contact')
        return F, Fdot, Fddot

    def reference(self, t):
        if self.base is not None and t < self.brake_start:
            return self.base.reference(t)
        plant, p, c = self.plant, self.p, self.cfg
        V, Q = self.plan.sol(t)
        D, A, pe, pg, y, _ = geom = plant.geometry(V)
        Fn, _, _ = self.load(t, V, Q, None, geom)
        qdot = (pe+pg+plant.head_Pa+Fn/A-plant.loss(Q)[0])/plant.I
        Fn, fdot, fddot = self.load(t, V, Q, qdot, geom)
        x = y+Fn/p['contact_stiffness']
        v = Q/A+fdot/p['contact_stiffness']
        acc = qdot/A+8*Q*Q/(math.pi*D**5)+fddot/p['contact_stiffness']
        # Positive-direction Coulomb compensation supplies a small starting
        # torque, not an invented initial contact preload. The forward plant
        # retains its original regularised friction law; it is not changed.
        friction = p['friction']
        torque = c.radius_m*(plant.M*acc+p['damping']*v+friction+Fn)
        if abs(torque) >= c.torque_cap:
            raise ValueError(f'Feedforward torque {torque:.6g} exceeds conditional envelope')
        theta = x/c.radius_m+math.asin(torque/c.torque_cap)/50
        count = theta*c.steps_per_rev/(2*math.pi)
        return dict(V=V,Q=Q,Fn=Fn,x=x,v=v,acc=acc,torque=torque,count=count)

    def scalar(self, t):
        if self.base is not None and t < self.brake_start:
            return self.base.scalar(t)
        return float(self.command(t)), float(self.command(t, 1))

    def acceleration(self, t):
        if self.base is not None and t < self.brake_start:
            return self.base.acceleration(t)
        return float(self.command(t, 2))

    def pulse_state(self, t):
        if np.ndim(t) == 0:
            return self.scalar(float(t))
        rows = np.array([self.scalar(float(tt)) for tt in t])
        return rows[:, 0], rows[:, 1]


class ClosedValveBrake:
    """Continue the motor phase smoothly after ideal valve closure.

    Rate and acceleration match at the switching instant. The brake does not
    prescribe actual rotor motion; the closed-valve mechanical ODE is solved.
    """
    enabled = True
    valve_open = False
    requested_steps_s = 0.0

    def __init__(self, before, start):
        self.before, self.start, self.cfg = before, start, before.cfg
        self.count0, self.rate0 = before.scalar(start)
        self.accel0 = before.acceleration(start)
        limit = self.cfg.closed_brake_accel
        if self.rate0 < 0 or abs(self.accel0) >= limit:
            raise ValueError('Invalid terminal-brake initial command state')
        self.duration = max(.1, 1.8*self.rate0/limit)
        for _ in range(30):
            u = np.linspace(0., 1., 1001)
            rate, accel = self._rate_accel(u)
            if min(rate) < -1e-7:
                raise ValueError('Terminal brake requires negative STEP rate')
            if max(abs(accel)) <= limit:
                break
            self.duration *= 1.08
        else:
            raise ValueError('Could not satisfy terminal command acceleration bound')
        self.end = start+self.duration

    def _rate_accel(self, u):
        r, a, T = self.rate0, self.accel0, self.duration
        rate = r*(2*u**3-3*u*u+1)+T*a*(u**3-2*u*u+u)
        accel = r*(6*u*u-6*u)/T+a*(3*u*u-4*u+1)
        return rate, accel

    def scalar(self, t):
        if t < self.start:
            return self.before.scalar(t)
        u = min((t-self.start)/self.duration, 1.)
        T, r, a = self.duration, self.rate0, self.accel0
        count = self.count0+T*r*(.5*u**4-u**3+u)
        count += T*T*a*(.25*u**4-2*u**3/3+.5*u*u)
        return count, float(self._rate_accel(u)[0])

    def acceleration(self, t):
        if t < self.start:
            return self.before.acceleration(t)
        if t >= self.end:
            return 0.0
        return float(self._rate_accel((t-self.start)/self.duration)[1])

    def pulse_state(self, t):
        if np.ndim(t) == 0:
            return self.scalar(float(t))
        rows = np.array([self.scalar(float(tt)) for tt in t])
        return rows[:, 0], rows[:, 1]
