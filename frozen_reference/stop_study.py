"""Conditional 15 mL/s motor-input / instantaneous binary-valve simulation.

Version 0.2: load-aware startup, taut initial cable, no initial preload,
instantaneous valve switching, and a distinct downstream residual-water tail.
Actual motion and flow
are solved by the forward plant, not substituted from the reference schedule.
The old model is retained byte-for-byte in source_snapshot for regression.
All physical calculations are in SI units. Python, not MATLAB, executes this.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "source_snapshot"))
import reference_model as ref
from motor_feedforward import LoadAwareProfile, ClosedValveBrake


@dataclass(frozen=True)
class Settings:
    q_target: float = 15e-6
    rise_s: float = 2.0
    slowdown_s: float = 3.0
    q_close_design: float = 3e-6
    steps_per_rev: int = 3200
    radius_m: float = 0.0025
    max_rate: float = 3000.0
    max_accel: float = 6000.0
    closed_brake_accel: float = 2000.0
    tail_volume_mL: float = 2.0
    tail_compensation_mL: float = 2.0
    downstream_length_m: float = 0.33
    torque_cap: float = 0.29658517619749203
    hold_after_seal_s: float = 0.5
    rtol: float = 2e-8
    max_step: float = 0.025

    def __post_init__(self):
        if min(self.q_target,self.rise_s,self.slowdown_s,self.q_close_design,
               self.max_rate,self.max_accel,self.closed_brake_accel,
               self.torque_cap,self.max_step) <= 0:
            raise ValueError('Positive timing, flow, and actuator scenario limits required')
        if self.q_close_design >= self.q_target:
            raise ValueError('The slowdown flow design must be below the plateau')
        if min(self.tail_volume_mL,self.tail_compensation_mL,self.downstream_length_m) < 0:
            raise ValueError('Tail inventory and compensation must be nonnegative')


class LegacyMotorProfile:
    """Previous geometric command, retained ONLY for open-ODE regression.

    It is not used by the new startup or terminal simulations.
    """
    enabled = True
    valve_open = True  # actual hydraulic closure is supplied separately

    def __init__(self, p, cfg, brake_start=None):
        self.p, self.cfg, self.brake_start = p, cfg, brake_start
        self.m_per_step = 2 * math.pi * cfg.radius_m / cfg.steps_per_rev
        self.brake_end = math.inf
        if brake_start is not None:
            self.count0, self.rate0 = self.before(brake_start)
            self.brake_duration = 1.5 * self.rate0 / cfg.max_accel
            self.brake_end = brake_start + self.brake_duration

    def before(self, t):
        c, p = self.cfg, self.p
        s = min(max(float(t), 0.0), c.rise_s)
        q = 0.5*c.q_target*(1-math.cos(math.pi*s/c.rise_s))
        void = 0.5*c.q_target*(s-c.rise_s/math.pi*math.sin(math.pi*s/c.rise_s))
        void += c.q_target*max(float(t)-c.rise_s, 0.0)
        V = p['V0']-void
        if V <= p['Vmin']:
            raise ValueError('Command geometry went beyond the declared spherical domain')
        D = (6*(V+p['Vsilicone'])/math.pi)**(1/3)
        return math.pi*(p['D0']-D)/self.m_per_step, (2*q/D**2)/self.m_per_step

    def scalar(self, t):
        if self.brake_start is None or t <= self.brake_start:
            return self.before(t)
        s = min((t-self.brake_start)/self.brake_duration, 1.0)
        count = self.count0+self.rate0*self.brake_duration*(s-s**3+0.5*s**4)
        rate = self.rate0*(1-3*s*s+2*s**3)
        return count, max(rate, 0.0)

    def pulse_state(self, t):
        if np.ndim(t) == 0:
            return self.scalar(float(t))
        a = np.array([self.scalar(float(tt)) for tt in np.asarray(t)])
        return a[:, 0], a[:, 1]

    @property
    def requested_steps_s(self):
        return 0.0  # not used by the new simulation; legacy interface only


class Plant:
    """Scalar implementation of the existing ODE, checked against snapshot."""
    def __init__(self, cfg=Settings()):
        self.cfg = cfg
        self.p = ref.load_plant()
        self.p['head'] = +0.1
        self.h = ref.load_reference(self.p)
        self.a = ref.load_actuator('illustrative_test_profile.json', allow_illustrative=True)
        self.a.update(torque_Nm=[cfg.torque_cap]*2,
                      torque_basis='Hypothetical flat envelope, not measured 12 V running torque',
                      purpose='Conditional 15 mL/s startup and termination study, not hardware validation',
                      microstep_basis='1/16 is the declared simulation subdivision',
                      max_step_rate_s=cfg.max_rate,acceleration_steps_s2=cfg.max_accel,
                      supply_voltage_V=12.0,current_assumed_A=0.67,
                      current_limit_measured=False,current_regulator_simulated=False,
                      torque_curve_matches_drive_conditions=False)
        if not self.h['nonlinear_pipe'] or self.h['additional_component_K_Pa_s2_m6']:
            raise ValueError('Scalar regression core supports the current unchanged hydraulic configuration only')
        self.M = ref.step4.total_equiv_mass(self.a)
        self.I = self.h['Ih_Pa_s2_m3']
        self.re_factor = self.p['rho']*self.p['pipe_diameter']/(self.p['water_viscosity']*self.h['pipe_area_m2'])
        self.pipe_C = self.h['pipe_length_m']/self.p['pipe_diameter']*self.p['rho']/(2*self.h['pipe_area_m2']**2)
        self.head_Pa = self.p['rho']*self.p['gravity']*self.p['head']
        self.downstream_water_mL=self.h['pipe_area_m2']*cfg.downstream_length_m*1e6
        self.pipe_water_mL=self.h['pipe_area_m2']*self.h['pipe_length_m']*1e6
        if cfg.tail_volume_mL > self.downstream_water_mL:
            raise ValueError('Tail exceeds known downstream tube inventory; add measured component inventory first')

    def geometry(self, V):
        p = self.p
        ri = (3*V/(4*math.pi))**(1/3)
        D = (6*(V+p['Vsilicone'])/math.pi)**(1/3)
        A = D*D/2
        ua, ub = p['Dref_inner']/(2*ri), p['Dref_outer']/D
        pe = -2*p['mu']*(ua+ua**4/4-ub-ub**4/4)
        pg = p['cage_mass']*p['gravity']*p['gravity_height_slope']/(math.pi*A)
        dp = (2*p['mu']/3)*(ua*(1+ua**3)/V-ub*(1+ub**3)/(V+p['Vsilicone']))-2*pg/(3*(V+p['Vsilicone']))
        return D, A, pe, pg, math.pi*(p['D0']-D), dp

    def loss(self, Q):
        p, h = self.p, self.h
        aq = abs(Q)
        re = self.re_factor*aq
        Rlam = h['R_pipe_laminar_Pa_s_m3']
        loss, slope = p['Rh']*Q, p['Rh']
        if re > h['Re_laminar']:
            span = h['Re_turbulent']-h['Re_laminar']
            s = min((re-h['Re_laminar'])/span, 1.0)
            w, dw = s*s*(3-2*s), 6*s*(1-s)/span
            arg = 6.9/re+(h['pipe_absolute_roughness_m']/(3.7*p['pipe_diameter']))**1.11
            den = -1.8*math.log10(arg)
            f = den**-2
            df = -2*den**-3*(1.8*6.9/(math.log(10)*arg*re**2))
            turb = self.pipe_C*f*aq*aq
            dturb = self.pipe_C*(2*f*aq+df*self.re_factor*aq*aq)
            loss += math.copysign(w*(turb-Rlam*aq), Q)
            slope += dw*self.re_factor*(turb-Rlam*aq)+w*(dturb-Rlam)
        return loss, slope

    def obs(self, t, z, profile, *, sealed=False):
        p, c = self.p, self.cfg
        x, v, V, Q = z[:4]
        D, A, pe, pg, y, dpdV = self.geometry(V)
        overlap = max(x-y, 0.0)
        Fn = p['contact_stiffness']*overlap
        count, rate = profile.scalar(t)
        phi = 50*(2*math.pi*count/c.steps_per_rev-x/c.radius_m)
        torque = c.torque_cap*math.sin(phi)
        Fm = torque/c.radius_m
        fr = p['friction']*math.tanh(v/p['friction_velocity'])
        acc = (Fm-p['damping']*v-fr-Fn)/self.M
        tension = self.a['nonmotor_equiv_mass_kg']*acc+p['damping']*v+fr+Fn
        drive = pe+pg+Fn/A+self.head_Pa
        base_loss, base_slope = self.loss(Q)
        gate = 0.0 if sealed else 1.0
        loss, slope = base_loss, base_slope
        qdot = 0.0 if sealed else (drive-loss)/self.I
        return dict(D=D,A=A,pe=pe,pg=pg,y=y,dpdV=dpdV,overlap=overlap,Fn=Fn,
                    phi=phi,torque=torque,Fm=Fm,fr=fr,acc=acc,tension=tension,
                    drive=drive,loss=loss,slope=slope,Qdot=qdot,gate=gate,
                    rate=rate,count=count,pressure=drive-self.head_Pa)

    def rhs(self, t, z, profile, **kw):
        o = self.obs(t,z,profile,**kw)
        v,Q = z[1],z[3]
        power = o['Fm']*v+self.head_Pa*Q-self.p['damping']*v*v-o['fr']*v-o['loss']*Q
        return np.array([v,o['acc'],-Q,o['Qdot'],Q,power])

    def jac(self, t,z,profile,**kw):
        o = self.obs(t,z,profile,**kw)
        p,c = self.p,self.cfg
        x,v,V,Q = z[:4]
        A=o['A']; k=p['contact_stiffness'] if o['overlap']>0 else 0.0
        fnV=k/A
        dFm_dx=-c.torque_cap*math.cos(o['phi'])*50/c.radius_m**2
        dfr=p['friction']/p['friction_velocity']*(1-math.tanh(v/p['friction_velocity'])**2)
        J=np.zeros((6,6))
        J[0,1]=1.
        J[1,:3]=[(dFm_dx-k)/self.M,(-p['damping']-dfr)/self.M,-fnV/self.M]
        J[2,3]=-1.
        if not kw.get('sealed',False):
            dP=o['dpdV']+fnV/A-o['Fn']*(2/(math.pi*o['D']))/A**2
            J[3,[0,2,3]]=[k/A/self.I,dP/self.I,-o['slope']/self.I]
        J[4,3]=1.
        J[5,[0,1,3]]=[v*dFm_dx,o['Fm']-2*p['damping']*v-o['fr']-v*dfr,
                       self.head_Pa-o['loss']-Q*o['slope']]
        return J

    def guard_values(self,t,z,profile,**kw):
        o=self.obs(t,z,profile,**kw)
        return [math.pi/2-abs(o['phi']),self.p['xmax']-z[0],z[2]-self.p['Vmin'],
                z[0]+1e-6,2*math.pi-abs(z[1]/self.cfg.radius_m),o['tension']+1e-5,
                z[3]+1e-13]

    def integrate(self,tspan,z,profile,*,target=None,**kw):
        names=['phase_boundary','50mm_limit','sphere_limit','extension_boundary',
               'speed_envelope_boundary','cable_tension_boundary','reverse_flow_boundary']
        n_guard=6 if kw.get('sealed',False) else 7
        events=[]
        for idx in range(n_guard):
            def evt(t,z,idx=idx): return self.guard_values(t,z,profile,**kw)[idx]
            evt.terminal=True; evt.direction=-1
            events.append(evt)
        if target is not None:
            def volume_event(t,z): return target*1e-6-(self.p['V0']-z[2])
            volume_event.terminal=True; volume_event.direction=-1
            events.append(volume_event)
            names.append('target_discharge_reached')
        sol=solve_ivp(lambda t,z:self.rhs(t,z,profile,**kw),tspan,z,method='Radau',
                      rtol=self.cfg.rtol,atol=[2e-11,2e-10,2e-13,2e-14,2e-13,2e-10],
                      jac=lambda t,z:self.jac(t,z,profile,**kw),max_step=self.cfg.max_step,
                      dense_output=True,events=events)
        if not sol.success: raise RuntimeError(sol.message)
        stops=[names[i] for i,v in enumerate(sol.t_events) if len(v)]
        sol.stop=stops[0] if stops else 'segment_complete'
        return sol


def regression_checks(plant):
    p,h=plant.p,plant.h
    profile=LegacyMotorProfile(p,plant.cfg)
    max_rhs=max_jac=max_loss=0.
    for V in [500e-6,350e-6,200e-6,100e-6]:
        D,A,pe,pg,y,dp=plant.geometry(V)
        for Q in [0.,2e-6,8e-6,15e-6]:
            l,s=plant.loss(Q)
            max_loss=max(max_loss,abs(l-float(ref.hydraulic_loss(Q,p,h)['dp_loss'])),
                         abs(s-float(ref.hydraulic_slope(Q,p,h)))/p['Rh'])
            # Choose matching command phase at a valid time; state need not be a trajectory.
            z=np.array([y+0.0001,0.001,V,Q,0.,0.])
            a=plant.rhs(8.,z,profile)
            b=ref.rhs(8.,z,p,h,actuator=plant.a,segment=profile)
            j1=plant.jac(8.,z,profile)
            j2=ref.jacobian(8.,z,p,h,actuator=plant.a,segment=profile)
            max_rhs=max(max_rhs,float(np.max(abs(a-b)/(1+abs(b)))))
            max_jac=max(max_jac,float(np.max(abs(j1-j2)/(1+abs(j2)))))
    assert max_rhs<1e-8 and max_jac<1e-7 and max_loss<1e-6
    return dict(rhs_relative_discrepancy=max_rhs,jac_relative_discrepancy=max_jac,
                loss_discrepancy=max_loss)


def baseline(plant):
    profile=LoadAwareProfile(plant)
    z=np.array([0.,0.,plant.p['V0'],0.,0.,0.])
    sol=plant.integrate((0.,31.),z,profile,target=400.)
    if sol.stop!='target_discharge_reached':
        raise RuntimeError('Base case stopped at '+sol.stop)
    sol.profile=profile
    return sol,profile


def terminal_run(plant,base,trigger_mL,*,target_mL=None):
    """Slow the motor smoothly, then close the ideal valve at target volume.

    Without target_mL this is a calibration trial closing at the end of the
    3 s load ramp. With target_mL the actual forward-state volume event closes
    the valve. The subsequent motor brake retains finite inertia and contact.
    """
    cfg,p=plant.cfg,plant.p
    tb=brentq(lambda t:1e6*(p['V0']-base.sol(t)[2])-trigger_mL,
              cfg.rise_s,base.t[-1],xtol=1e-10)
    open_profile=LoadAwareProfile(plant,base=base.profile,brake_start=tb)
    end=tb+cfg.slowdown_s+(3. if target_mL is not None else 0.)
    valve_target=None if target_mL is None else target_mL-cfg.tail_compensation_mL
    open_sol=plant.integrate((tb,end),base.sol(tb),open_profile,target=valve_target)
    expected='target_discharge_reached' if target_mL is not None else 'segment_complete'
    if open_sol.stop != expected:
        raise RuntimeError('Open-valve slowdown stopped at '+open_sol.stop)
    tc=float(open_sol.t[-1])
    z=open_sol.y[:,-1].copy()
    q_left=float(z[3])
    tail_I=plant.I*cfg.downstream_length_m/plant.h['pipe_length_m'] if cfg.tail_volume_mL else 0.
    removed_energy=.5*(plant.I-tail_I)*q_left*q_left
    if cfg.tail_volume_mL and q_left <= 0:
        raise ValueError('A positive assumed tail requires positive flow at closure')
    tail_duration=2*cfg.tail_volume_mL*1e-6/q_left if cfg.tail_volume_mL else 0.
    volume_before_reset=float(z[2])
    # User-selected ideal hybrid reset THROUGH THE VALVE. It neither emits
    # nor removes water at the switching instant. The downstream water keeps
    # its initial velocity and subsequently drains through a separate, finite
    # inventory tail. Only upstream kinetic energy is removed at this reset.
    # The associated pressure impulse / water hammer is NOT resolved here.
    z[3]=0.; z[5]-=removed_energy
    profile=ClosedValveBrake(open_profile,tc)
    tm=profile.end
    outlet_stop=tc+tail_duration
    closed_sol=plant.integrate((tc,max(tm,outlet_stop)+cfg.hold_after_seal_s),z,profile,sealed=True)
    if closed_sol.stop != 'segment_complete':
        raise RuntimeError('Closed-valve motor brake stopped at '+closed_sol.stop)
    pieces=[open_sol,closed_sol]
    Vfinal=closed_sol.y[2,-1]
    return dict(stop='segment_complete',trigger_mL=trigger_mL,brake_start_s=tb,
                close_command_s=tc,motor_stop_command_s=tm,seal_s=tc,
                cycle_end_s=max(tm,outlet_stop),valve_delay_s=0.,valve_travel_s=0.,
                final_mL=1e6*(p['V0']-Vfinal)+cfg.tail_volume_mL,
                bladder_discharge_mL=1e6*(p['V0']-Vfinal),
                volume_at_close_command_mL=1e6*(p['V0']-volume_before_reset),
                post_close_command_mL=1e6*(volume_before_reset-Vfinal),
                outlet_tail_mL=cfg.tail_volume_mL,
                tail_compensation_mL=cfg.tail_compensation_mL,
                tail_duration_s=tail_duration,outlet_stop_s=outlet_stop,
                downstream_inertance=tail_I,
                flow_before_close_mL_s=q_left*1e6,
                closure_removed_fluid_energy_J=removed_energy,
                close_by_volume_event=target_mL is not None,
                pieces=pieces,profile=profile)


def tail_response(t,case):
    """Finite-volume phenomenological tail, not a calibrated pipe-drain law.

    It starts at the pre-close outlet flow and ends smoothly after T, with
    integral equal to the assumed residual volume: T=2*V_tail/Q(close-).
    No valve leakage is allowed; all tail water comes from pipe inventory.
    """
    t=np.asarray(t,dtype=float)
    T=case['tail_duration_s']
    if T == 0:
        return np.zeros_like(t),np.zeros_like(t)
    elapsed=np.clip(t-case['seal_s'],0.,T)
    q0=case['flow_before_close_mL_s']*1e-6
    q=.5*q0*(1+np.cos(np.pi*elapsed/T))
    drained=.5*q0*(elapsed+T/np.pi*np.sin(np.pi*elapsed/T))
    q=np.where(t>=case['seal_s'],q,0.)
    return q,drained


def trajectory(plant,base,case):
    tb=case['brake_start_s']
    ts=np.unique(np.r_[np.arange(0.,tb,.025),base.t[base.t<tb],
                        base.profile.grid[base.profile.grid<tb],tb])
    arrays=[base.sol(ts)]; time_arrays=[ts]; mode_arrays=[np.ones(len(ts))]
    for idx,piece in enumerate(case['pieces']):
        command_grid=case['profile'].before.grid if idx==0 else np.array([])
        command_grid=command_grid[(command_grid>=piece.t[0])&(command_grid<=piece.t[-1])]
        t=np.unique(np.r_[piece.t,command_grid,np.arange(piece.t[0],piece.t[-1],.025),piece.t[-1]])
        time_arrays.append(t); arrays.append(piece.sol(t)); mode_arrays.append(np.full(len(t),1.-idx))
    t=np.concatenate(time_arrays); z=np.concatenate(arrays,axis=1)
    modes=np.concatenate(mode_arrays)
    # Keep left/right limits at ideal closure to plot the vertical Q jump.
    keep=np.r_[(np.diff(t)>1e-12)|(np.diff(modes)!=0),True]
    t,z,modes=t[keep],z[:,keep],modes[keep]
    obs=[plant.obs(float(tt),zz,case['profile'],sealed=mode==0)
         for tt,zz,mode in zip(t,z.T,modes)]
    out={k:np.array([o[k] for o in obs]) for k in obs[0]}
    qt,drained=tail_response(t,case)
    closed=modes==0
    qt=np.where(closed,qt,0.)
    tail_energy=.5*case['downstream_inertance']*qt**2
    tail_E0=.5*case['downstream_inertance']*(case['flow_before_close_mL_s']*1e-6)**2
    tail_dissipation=np.where(closed,tail_E0-tail_energy,0.)
    E=(0.5*plant.M*z[1]**2+0.5*plant.p['contact_stiffness']*out['overlap']**2
       +ref.step4.legacy_baseline.shell_energy(z[2],plant.p)
       +plant.p['cage_mass']*plant.p['gravity']*plant.p['gravity_height_slope']*(out['D']-plant.p['D0'])
       +0.5*plant.I*z[3]**2+tail_energy)
    out.update(t=t,z=z,energy_J=E,voided_mL=1e6*(plant.p['V0']-z[2]),
               outlet_flow_m3_s=np.where(closed,qt,z[3]),
               outlet_collected_mL=1e6*(plant.p['V0']-z[2]+drained),
               tail_drained_mL=drained*1e6,tube_water_mL=plant.pipe_water_mL-drained*1e6,
               tail_energy_J=tail_energy,tail_dissipation_J=tail_dissipation,
               energy_work_J=z[5]-tail_dissipation,
               command_accel=np.array([case['profile'].acceleration(tt) for tt in t]))
    return out


def diagnose(plant,data):
    z,t=data['z'],data['t']
    cfg=plant.cfg
    mass_error=float(np.max(abs(z[4]-(plant.p['V0']-z[2])))*1e6)
    energy_error=float(np.max(abs(data['energy_J']-data['energy_J'][0]-data['energy_work_J'])))
    total_mass_error=float(np.max(abs(plant.p['V0']*1e6+plant.pipe_water_mL-
                                      z[2]*1e6-data['tube_water_mL']-data['outlet_collected_mL'])))
    result=dict(max_flow_mL_s=float(max(z[3])*1e6),max_torque_Nm=float(max(abs(data['torque']))),
                max_motor_rpm=float(max(abs(z[1]))/cfg.radius_m*60/(2*math.pi)),
                max_command_steps_s=float(max(data['rate'])),
                max_command_acceleration_steps_s2=float(max(abs(data['command_accel']))),
                max_phase_deg=float(max(abs(data['phi']))*180/math.pi),
                minimum_cable_tension_N=float(min(data['tension'])),
                max_radial_contact_gap_um=float(max(np.maximum(data['y']-z[0],0.))/(2*math.pi)*1e6),
                minimum_contact_force_after_0p1s_N=float(min(data['Fn'][t>=.1])),
                final_D_mm=float(data['D'][-1]*1000),
                mass_balance_error_mL=mass_error,total_inventory_error_mL=total_mass_error,
                energy_balance_error_J=energy_error)
    assert mass_error<1e-6 and energy_error<1e-5
    assert total_mass_error<1e-6
    assert result['max_command_steps_s']<=cfg.max_rate+1e-6
    assert result['max_command_acceleration_steps_s2']<=cfg.max_accel+1e-3
    return result


def clean(case):
    return {k:v for k,v in case.items() if k not in ('pieces','profile')}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',action='store_true')
    args=parser.parse_args()
    plant=Plant()
    print('REGRESSION',regression_checks(plant),flush=True)
    tick=time.monotonic()
    base,profile=baseline(plant)
    print('BASE',base.stop,base.t[-1],'wall_s',time.monotonic()-tick,flush=True)
    if args.probe:
        case=terminal_run(plant,base,373.)
        print('PROBE',json.dumps(clean(case)),flush=True)
        print('DIAGNOSTICS',diagnose(plant,trajectory(plant,base,case)),flush=True)
        return
    results=ROOT/'results'
    results.mkdir(exist_ok=True)
    reports=[]; saved=[]
    for target,guess in [(200.,171.),(300.,271.),(400.,371.)]:
        # Offline nominal feedforward calibration, not a fit to experimental
        # data. Trigger on actual cumulative volume in subsequent scenarios.
        trigger=guess
        trials=[]
        for iteration in range(5):
            case=terminal_run(plant,base,trigger)
            if case['stop']!='segment_complete':
                raise RuntimeError(json.dumps(clean(case)))
            error=case['final_mL']-target
            trials.append(dict(trigger_mL=trigger,final_mL=case['final_mL']))
            print('NOMINAL',target,iteration,json.dumps(clean(case)),flush=True)
            if abs(error)<0.002: break
            # Braking travel changes mildly with trigger volume. A secant
            # update uses simulation, never an invented experimental label.
            slope=1.0
            if len(trials)>1:
                prev=trials[-2]
                slope=(case['final_mL']-prev['final_mL'])/(trigger-prev['trigger_mL'])
            trigger-=error/slope
        if abs(error)>=0.002: raise RuntimeError('Nominal trigger did not converge')
        case=terminal_run(plant,base,trigger,target_mL=target)
        error=case['final_mL']-target
        data=trajectory(plant,base,case)
        report=clean(case)
        report.update(target_mL=target,endpoint_error_mL=error,diagnostics=diagnose(plant,data),
                      nominal_trigger_calibration=trials)
        mask=(data['t']>=plant.cfg.rise_s+.05)&(data['t']<case['brake_start_s'])
        report['plateau_flow_min_max_mL_s']=[float(min(data['z'][3,mask])*1e6),
                                            float(max(data['z'][3,mask])*1e6)]
        np.savez_compressed(results/f'target_{int(target)}mL.npz',**data)
        cols=np.column_stack([data['t'],data['z'][0],data['z'][1],data['z'][2],data['z'][3],
                              data['voided_mL'],data['torque'],data['rate'],data['gate'],
                              data['Fn'],data['tension'],data['pressure'],data['phi'],
                              data['count']*2*math.pi/plant.cfg.steps_per_rev,data['command_accel'],
                              data['outlet_flow_m3_s'],data['outlet_collected_mL'],data['tube_water_mL']])
        np.savetxt(results/f'target_{int(target)}mL.csv',cols,delimiter=',',
                   header='time_s,cable_x_m,cable_v_m_s,water_V_m3,valve_Q_m3_s,bladder_discharge_mL,torque_Nm,command_steps_s,valve_open,contact_force_N,cable_tension_N,bladder_pressure_Pa,electrical_phase_rad,command_angle_rad,command_accel_steps_s2,outlet_Q_m3_s,outlet_collected_mL,tube_water_mL',comments='')
        reports.append(report); saved.append(data)
        print('VERIFIED_CASE',json.dumps(report),flush=True)
    startup=[]
    for tt in [0.,.01,.05,.1,.5,1.,1.5,2.,3.]:
        zz=base.sol(tt);o=plant.obs(tt,zz,profile)
        startup.append(dict(time_s=tt,flow_mL_s=zz[3]*1e6,contact_force_N=o['Fn'],
                            cable_tension_N=o['tension'],torque_Nm=o['torque'],
                            actual_motor_rpm=zz[1]/plant.cfg.radius_m*60/(2*math.pi)))
    summary=dict(model_version='0.2',settings=asdict(plant.cfg),plant=plant.p,hydraulics=plant.h,
                 actuator=plant.a,regression=regression_checks(plant),cases=reports,
                 startup_samples=startup,
                 tail_basis='Approximately 2 mL reported by user for this experiment; duration/shape unmeasured',
                 tail_inventory_known_downstream_mL=plant.downstream_water_mL,
                 tail_sensitivity=[dict(actual_tail_mL=v,final_volume_error_mL=v-plant.cfg.tail_compensation_mL)
                                   for v in [1.,1.5,2.,2.5,3.]],
                 data_kind='conditional_simulation_not_experiment',
                 hardware_validation=False,independent_experimental_validation=False,
                 assumptions=['15 mL/s is a design reference, actual Q is an ODE state',
                              'motor and full-open valve start together at t=0; cable initially taut',
                              'initial x=v=Q=0, V=500 mL, initial contact preload=0 N',
                              'load-aware inverse feedforward plus separate forward-dynamics replay',
                              'startup load ramp 2 s; bounded average microstep phase, not discrete pulses',
                              '6000 steps/s2 command limit is a revised design prior, not a measured limit',
                              'flat 0.296585 N m envelope is hypothetical, not a 12 V measured curve',
                              'outlet collection includes the approximately 2 mL downstream tail reported by the user',
                              'instantaneous closure at target minus estimated tail; valve flow reset to zero',
                              'tail duration follows 2*V_tail/Q(close-); tail shape is phenomenological',
                              'known downstream pipe inventory is reduced by the tail, not bladder volume',
                              'ideal closure removes upstream kinetic energy but no water volume',
                              'no distributed water hammer or tube-wall compliance',
                              'spherical shell assumed, including below natural volume'])
    (results/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    plot_summary(saved,reports,results)
    plot_startup(saved[-1],results)
    plot_closure(saved[-1],reports[-1],results)
    print('DONE wall_s',time.monotonic()-tick,flush=True)


def plot_summary(data,reports,folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold'})
    fig,axs=plt.subplots(2,1,figsize=(11,7.8),sharex=True)
    fig.subplots_adjust(left=.095,right=.985,top=.845,bottom=.15,hspace=.15)
    colors=['#1c7f85','#e69031','#596bb2']
    for d,r,color in zip(data,reports,colors):
        label=f"{int(r['target_mL'])} mL target"
        axs[0].plot(d['t'],d['outlet_flow_m3_s']*1e6,color=color,lw=2,label=label)
        axs[1].plot(d['t'],d['outlet_collected_mL'],color=color,lw=2)
        axs[0].plot(r['outlet_stop_s'],0.,'o',color=color,ms=5)
        axs[1].plot(r['outlet_stop_s'],r['final_mL'],'o',color=color,ms=5)
        axs[1].axhline(r['target_mL'],color=color,lw=.7,ls=':',alpha=.5)
    axs[0].axhline(15,color='#889099',ls='--',lw=1,label='15 mL/s design plateau')
    axs[0].set(ylabel='Final outlet flow (mL/s)',ylim=(-.5,17))
    axs[1].set(ylabel='Collected outlet volume (mL)',xlabel='Time (s)',ylim=(0,430))
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.54,.905),
               ncol=4,frameon=False,fontsize=9)
    for ax in axs: ax.grid(axis='y',alpha=.18)
    fig.suptitle('Immediate loading + ideal valve + residual-water tail',fontsize=16,y=.977)
    fig.text(.54,.935,'ALBA | 500 mL initial volume | 15 mL/s design plateau | +10 cm outlet head',
             ha='center',fontsize=11,color='#505762')
    tail=reports[0]['outlet_tail_mL']
    fig.text(.095,.03,f'Motor and valve start at t = 0. Valve closes at target minus {tail:g} mL; the outlet then drains the residual water.\nTail volume ~{tail:g} mL reported by user; duration unmeasured. Conditional 0.297 N m torque envelope, not hardware validation.',fontsize=9,color='#505762')
    fig.savefig(folder/'ALBA_15mLps_full_cycle.png',dpi=170)
    fig.savefig(folder/'ALBA_15mLps_full_cycle.svg')
    plt.close(fig)


def plot_startup(d,folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(3,1,figsize=(10,8),sharex=True)
    fig.subplots_adjust(left=.115,right=.97,top=.88,bottom=.14,hspace=.23)
    mask=d['t']<=3.5
    t=d['t'][mask]
    axs[0].plot(t,d['z'][3,mask]*1e6,color='#596bb2',lw=2)
    axs[0].axhline(15,color='#999999',ls='--',lw=1)
    axs[0].set(ylabel='Flow (mL/s)',ylim=(-.3,16))
    axs[1].plot(t,d['Fn'][mask],color='#1c7f85',lw=2)
    axs[1].set(ylabel='Contact load (N)',ylim=(-2,120))
    axs[2].plot(t,d['z'][1,mask]/.0025*60/(2*np.pi),color='#e69031',lw=2)
    axs[2].set(ylabel='Actual motor (rpm)',xlabel='Time after simultaneous valve / motor start (s)')
    for ax in axs:
        ax.axvline(2.,color='#b6bbc1',ls=':',lw=1)
        ax.grid(axis='y',alpha=.18)
    fig.suptitle('Startup detail: no seconds-long unloaded interval',fontsize=15,y=.965)
    fig.text(.54,.917,'Instant full opening | No initial slack or added preload | 2 s active load ramp',
             ha='center',color='#505762',fontsize=10)
    fig.text(.115,.035,'The plotted flow, contact force and rotor speed are forward-ODE outputs, not prescribed output curves.\nWater inertance makes Q start from zero; contact loading begins immediately for t > 0.',fontsize=9,color='#505762')
    fig.savefig(folder/'ALBA_15mLps_startup_detail.png',dpi=170)
    fig.savefig(folder/'ALBA_15mLps_startup_detail.svg')
    plt.close(fig)


def plot_closure(d,r,folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,1,figsize=(10,6.6),sharex=True)
    fig.subplots_adjust(left=.12,right=.97,top=.85,bottom=.17,hspace=.3)
    dt=d['t']-r['seal_s']
    mask=(dt>=-.25)&(dt<=r['tail_duration_s']+.25)
    axs[0].plot(dt[mask],d['outlet_flow_m3_s'][mask]*1e6,color='#596bb2',lw=2.3,label='Final outlet')
    axs[0].plot(dt[mask],d['z'][3,mask]*1e6,color='#e69031',lw=1.8,ls='--',label='Through valve')
    axs[0].set(ylabel='Flow (mL/s)',ylim=(-.15,3.5))
    axs[0].legend(loc='upper right',frameon=False)
    axs[1].plot(dt[mask],d['outlet_collected_mL'][mask],color='#596bb2',lw=2.3,label='Collected at outlet')
    axs[1].plot(dt[mask],d['voided_mL'][mask],color='#e69031',lw=1.8,ls='--',label='Displaced from bladder')
    axs[1].set(ylabel='Volume (mL)',xlabel='Time relative to valve closure (s)',
               ylim=(r['target_mL']-r['outlet_tail_mL']-.9,r['target_mL']+.2))
    axs[1].legend(loc='lower right',frameon=False,fontsize=9)
    for ax in axs:
        ax.axvline(0.,color='#68717c',ls=':',lw=1)
        ax.grid(axis='y',alpha=.18)
    axs[0].text(.025,3.32,'Valve fully closed',fontsize=9,color='#505762')
    fig.suptitle('Valve flow stops immediately; outlet water drains afterward',fontsize=14,y=.965)
    fig.text(.54,.907,f"400 mL target | Residual volume ~{r['outlet_tail_mL']:g} mL reported by user; duration unmeasured",
             ha='center',fontsize=11,color='#505762')
    fig.text(.12,.03,'Bladder water stays fixed after closure. The extra outlet water comes only from the downstream tube.\nResidual volume uses this experiment; the smooth tail duration and shape remain modelling assumptions.',fontsize=9,color='#505762')
    fig.savefig(folder/'ALBA_15mLps_closure_detail.png',dpi=170)
    fig.savefig(folder/'ALBA_15mLps_closure_detail.svg')
    plt.close(fig)


if __name__=='__main__':
    main()
