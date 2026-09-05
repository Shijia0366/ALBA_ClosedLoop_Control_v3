"""ALBA force-input reduced ODE. SI units throughout.

Physical states: x (cage/cable retraction), v, V (water volume).
Q and outlet-referenced pressure are algebraic outputs.  The spherical
silicone shell has its own volume-dependent energy; unilateral contact
couples it to the cage without forcing a shrinking bladder to pull the cage.

This is a provisional model, not a calibrated reproduction of Fixed-ID.
The MATLAB implementation is in ../matlab. No ideal velocity source is used.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad, solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reference_results"


def parameters():
    p = dict(V0=500e-6, Dref_outer=.070, Dref_inner=.066,
             Dmin=.050, cable_gap=.003, material="Ecoflex 00-20",
             material_confirmed_by_user=True,
             stress100=8*6894.757293168, rho=1000., gravity=9.81,
             water_viscosity=1.0e-3, pipe_diameter=.004,
             pipe_segments=[.03,.03,.30], head=-.10,
             shaft_diameter=.005, effective_radius=.0025,
             mass_equivalent=.50, damping=20., friction=.10,
             friction_velocity=1e-5, force_max=15., contact_stiffness=2e5,
             Rh=4954212868.228043, force_ramp=1., time_horizon=3600.,
             rtol=2e-8, atol=[2e-11,2e-10,2e-13,2e-13,2e-11],
             max_step=5., pressure_reference="lumped bladder outlet gauge pressure")
    p["mu"] = p["stress100"]/(2.-2.**-2)
    p["Vsilicone"] = np.pi/6*(p["Dref_outer"]**3-p["Dref_inner"]**3)
    p["Vnatural"] = np.pi/6*p["Dref_inner"]**3
    p["D0"] = (6*(p["V0"]+p["Vsilicone"])/np.pi)**(1/3)
    p["L0"] = np.pi*(p["D0"]+2*p["cable_gap"])
    p["xmax"] = np.pi*(p["D0"]-p["Dmin"])
    p["Vmin"] = np.pi/6*p["Dmin"]**3-p["Vsilicone"]
    p["pipe_length"] = sum(p["pipe_segments"])
    p["pipe_area"] = np.pi*p["pipe_diameter"]**2/4
    p["R_tube"] = 128*p["water_viscosity"]*p["pipe_length"]/(np.pi*p["pipe_diameter"]**4)
    p["I_tube"] = p["rho"]*p["pipe_length"]/p["pipe_area"]
    p["tube_volume"] = p["pipe_area"]*p["pipe_length"]
    return p


def shell(V, p):
    """Uniform spherical incompressible Neo-Hookean shell, analytic pressure.

    pe = dU/dV.  Pressure is an effective lumped value assigned to the outlet;
    internal vertical pressure gradients and the bladder-neck geometry are
    omitted. Below the natural volume this is a constrained-sphere
    compression extrapolation, not a buckling/collapse law.
    """
    V = np.asarray(V, dtype=float)
    if np.any(V <= 0) or not np.all(np.isfinite(V)):
        raise ValueError("Positive finite water volume is required")
    a = (3*V/(4*np.pi))**(1/3)
    b = (3*(V+p["Vsilicone"])/(4*np.pi))**(1/3)
    a0, b0 = p["Dref_inner"]/2, p["Dref_outer"]/2
    ua, ub = a0/a, b0/b
    pe = -2*p["mu"]*(ua+ua**4/4-ub-ub**4/4)
    D = 2*b
    area = D*D/2  # -dV/dy, y=pi*(D0-D)
    return dict(D=D, Di=2*a, wall=b-a, A=area, pe=pe,
                y=np.pi*(p["D0"]-D), Fel=-area*pe,
                below_natural=V < p["Vnatural"])


def shell_energy(V, p):
    """Independent quadrature used only for energy verification/output."""
    V = np.asarray(V, dtype=float)
    a0, b0 = p["Dref_inner"]/2, p["Dref_outer"]/2
    z,w = np.polynomial.legendre.leggauss(48)
    R = (a0+b0)/2+(b0-a0)/2*z
    w = w*(b0-a0)/2
    a3 = 3*V/(4*np.pi)
    r = (R**3-a0**3+a3[...,None])**(1/3)
    W = p["mu"]/2*((R/r)**4+2*(r/R)**2-3)
    return 4*np.pi*np.sum(w*W*R**2,axis=-1)


def force_input(t, requested_force, p):
    """Bounded FORCE input, not step-frequency or prescribed motion."""
    level = np.clip(requested_force, 0., p["force_max"])
    t = np.asarray(t, float)
    if p["force_ramp"] == 0:
        return level*np.ones_like(t)
    s = np.clip(t/p["force_ramp"],0.,1.)
    return level*(1-np.cos(np.pi*s))/2


def outputs(t, state, requested_force, p):
    x,v,V = np.asarray(state)[:3]
    a = shell(V,p)
    delta = np.maximum(x-a["y"],0.)
    contact = p["contact_stiffness"]*delta
    pressure = a["pe"]+contact/a["A"]
    drive = pressure+p["rho"]*p["gravity"]*p["head"]
    # Discharge-only approximation: no external reservoir supplies reverse
    # liquid. Does not model return of water stored in the pipe, air ingestion,
    # menisci, rotor breakaway or de-priming.
    Q = np.maximum(drive,0.)/p["Rh"]
    friction = p["friction"]*np.tanh(v/p["friction_velocity"])
    Fm = force_input(t,requested_force,p)
    a.update(x=x,v=v,V=V,voided=p["V0"]-V,Fn=contact,p=pressure,
             drive=drive,Q=Q,Fm=Fm,Ffriction=friction,
             overlap=delta,radial_gap=(a["y"]-x)/(2*np.pi),
             Dcage=p["D0"]-x/np.pi,Lring=p["L0"]-x,
             torque=Fm*p["effective_radius"],
             ydot=Q/a["A"],Re=p["rho"]*Q*p["pipe_diameter"]/(p["water_viscosity"]*p["pipe_area"]))
    return a


def rhs(t, state, requested_force, p):
    """Three physical states; optional two independent audit integrals."""
    a = outputs(t,state,requested_force,p)
    acceleration = (a["Fm"]-p["damping"]*a["v"]-a["Ffriction"]-a["Fn"])/p["mass_equivalent"]
    physical = [a["v"], acceleration, -a["Q"]]
    if len(state) == 3:
        return physical
    power = (a["Fm"]*a["v"]+p["rho"]*p["gravity"]*p["head"]*a["Q"]
             -p["damping"]*a["v"]**2-a["Ffriction"]*a["v"]-p["Rh"]*a["Q"]**2)
    return physical+[a["Q"],power]


def jacobian(t,state,requested_force,p):
    """Analytic stiff-system Jacobian; avoids perturbing volume across contact."""
    a=outputs(t,state,requested_force,p)
    V=float(state[2]); v=float(state[1])
    ri=(3*V/(4*np.pi))**(1/3)
    ro=(3*(V+p["Vsilicone"])/(4*np.pi))**(1/3)
    ua=p["Dref_inner"]/(2*ri); ub=p["Dref_outer"]/(2*ro)
    dp=(2*p["mu"]/3)*(ua*(1+ua**3)/V-ub*(1+ub**3)/(V+p["Vsilicone"]))
    k=p["contact_stiffness"] if a["overlap"]>0 else 0.
    area=float(a["A"]); darea=2/(np.pi*a["D"])
    fnx=k; fnV=k/area
    qx=qV=0.
    if a["drive"]>0:
        qx=k/(area*p["Rh"])
        qV=(dp+fnV/area-a["Fn"]*darea/area**2)/p["Rh"]
    df=p["friction"]/p["friction_velocity"]*(1-np.tanh(v/p["friction_velocity"])**2)
    J=np.zeros((len(state),len(state)))
    J[0,1]=1
    J[1,:3]=[-fnx/p["mass_equivalent"],-(p["damping"]+df)/p["mass_equivalent"],-fnV/p["mass_equivalent"]]
    J[2,:3]=[-qx,0,-qV]
    if len(state)==5:
        J[3,:3]=[qx,0,qV]
        hydraulic_gradient=p["rho"]*p["gravity"]*p["head"]-2*p["Rh"]*a["Q"]
        J[4,:3]=[hydraulic_gradient*qx,a["Fm"]-a["Ffriction"]-v*df-2*p["damping"]*v,hydraulic_gradient*qV]
    return J


def run_case(name, force, p, horizon=None, initial=None):
    p = dict(p)
    end = p["time_horizon"] if horizon is None else float(horizon)
    if initial is None:
        initial = [0.,0.,p["V0"]]
    if len(initial) != 3 or initial[2] <= p["Vmin"] or initial[0] >= p["xmax"]:
        raise ValueError("Initial state must be strictly within the model domain")
    def cable_limit(t,z): return p["xmax"]-z[0]
    def spherical_domain(t,z): return z[2]-p["Vmin"]
    def extension_guard(t,z): return z[0]+1e-6
    def free_shape_limit(t,z):
        # A detached sphere below its natural volume needs a collapse law.
        # Stop before using a smooth compressed free sphere to predict a
        # spurious passive-drainage equilibrium. The force tolerance is only
        # for domain classification, not a fitted contact constitutive law.
        fn=p["contact_stiffness"]*max(z[0]-float(shell(z[2],p)["y"]),0.)
        return max(z[2]-p["Vnatural"],1e-6*(fn-1e-5))
    for event in (cable_limit,spherical_domain,extension_guard,free_shape_limit):
        event.terminal=True
        event.direction=-1
    sol=solve_ivp(lambda t,z:rhs(t,z,force,p),(0.,end),list(initial)+[0.,0.],
                  method="Radau",rtol=p["rtol"],atol=p["atol"],
                  max_step=p["max_step"],dense_output=True,
                  jac=lambda t,z:jacobian(t,z,force,p),
                  events=(cable_limit,spherical_domain,extension_guard,free_shape_limit))
    if not sol.success:
        raise RuntimeError(sol.message)
    tend=sol.t[-1]
    # Include dense early samples to integrate short force/contact transients.
    t=np.unique(np.r_[np.linspace(0,min(3.,tend),1501),
                      np.arange(3.,tend,.25),sol.t,tend])
    z=sol.sol(t)
    a=outputs(t,z,force,p)
    a.update(t=t,z=z,energy=.5*p["mass_equivalent"]*z[1]**2+shell_energy(z[2],p)
             +.5*p["contact_stiffness"]*a["overlap"]**2)
    stop="time_horizon"
    for event,label in zip(sol.t_events,("cage_50mm_stop","spherical_50mm_domain_boundary","extension_model_boundary","free_collapse_model_required")):
        if len(event): stop=label
    targets={}
    for volume in (200.,300.,400.):
        Vtarget=p["V0"]-volume*1e-6
        if z[2,0]<=Vtarget:
            targets[str(int(volume))]=0.
        elif z[2,-1]<=Vtarget:
            targets[str(int(volume))]=float(brentq(lambda tt:sol.sol(tt)[2]-Vtarget,0.,tend,xtol=1e-9))
        else:
            targets[str(int(volume))]=None
    report=dict(name=name,requested_force_N=force,force_cap_N=p["force_max"],
                head_m=p["head"],Rh_Pa_s_m3=p["Rh"],duration_s=float(tend),stop=stop,
                voided_mL=float(a["voided"][-1]*1e6),remaining_mL=float(z[2,-1]*1e6),
                x_mm=float(z[0,-1]*1e3),final_outer_diameter_mm=float(a["D"][-1]*1e3),
                peak_flow_mL_s=float(np.max(a["Q"])*1e6),
                final_flow_mL_s=float(a["Q"][-1]*1e6),
                terminal_value_is_physical_final_emptying=False,
                final_contact_force_N=float(a["Fn"][-1]),
                final_radial_gap_mm=float(a["radial_gap"][-1]*1e3),
                max_contact_penetration_mm=float(np.max(a["overlap"])/(2*np.pi)*1e3),
                max_reynolds=float(np.max(a["Re"])),
                target_times_s=targets,
                entered_unvalidated_compression_branch=bool(np.any(a["below_natural"])))
    return dict(p=p,sol=sol,data=a,summary=report)


def estimate_resistance(p):
    """One conditional scalar fit, excluding generated R3 and terminal drip.

    With Fm=0 the bladder detaches and x remains zero. For a taut liquid
    column and the chosen shell prior, t(V)=Rh*integral(dV/(pe+rho*g*h)).
    Equal weight for each measured time point. This estimates a lumped
    dissipation/time scale, not separate valve/flowmeter/pipe coefficients.
    """
    src=json.loads((ROOT/"input_data/fixed_id_zero_load_excerpt.json").read_text())
    candidates=[]
    records=[]
    for trial in src["trials"]:
        S=[]
        for out in trial["volume_mL"]:
            lower=p["V0"]-out*1e-6
            value=quad(lambda V:1/(float(shell(V,p)["pe"])+p["rho"]*p["gravity"]*p["head"]),
                       lower,p["V0"],epsabs=1e-14,epsrel=1e-10)[0]
            S.append(value)
        S=np.asarray(S); observed=np.asarray(trial["elapsed_s"])
        R=float(S@observed/(S@S))
        candidates.append((S,observed))
        records.append(dict(trial=trial["id"],conditional_Rh=R))
    ss=np.concatenate([s for s,t in candidates]); tt=np.concatenate([t for s,t in candidates])
    pooled=float(ss@tt/(ss@ss))
    return dict(Rh=pooled,per_trial=records,fit_points=len(tt),
                time_residual_RMSE_s=float(np.sqrt(np.mean((pooled*ss-tt)**2))),
                source=src["source_name"],source_version=src["source_version"],
                scope="Conditional on shell pressure prior and actual -0.10m head. No independent pressure measurements; not uniquely identified plumbing resistance.",
                data_range="0g R1/R2, 25-250mL; R3 and terminal drip excluded")


def main():
    OUT.mkdir(exist_ok=True)
    p=parameters()
    fit=estimate_resistance(p)
    p["Rh"]=fit["Rh"]
    specs=[("raised_zero",0.,-.10), ("raised_100g_force",.981,-.10),
           ("raised_200g_force",1.962,-.10), ("lowered_zero",0.,.10),
           ("raised_force_cap",15.,-.10)]
    cases={}
    for name,F,head in specs:
        cases[name]=run_case(name,F,dict(p,head=head))
        print(json.dumps(cases[name]["summary"],ensure_ascii=False),flush=True)
    from verify_step3 import verify_all
    verification=verify_all(p,cases)
    (OUT/"parameters_SI.json").write_text(json.dumps(p,indent=2,ensure_ascii=False))
    (OUT/"resistance_estimate.json").write_text(json.dumps(fit,indent=2,ensure_ascii=False))
    (OUT/"summary.json").write_text(json.dumps({k:c["summary"] for k,c in cases.items()},indent=2))
    (OUT/"verification.json").write_text(json.dumps(verification,indent=2))
    from scipy.io import savemat
    matcases={}
    for name,c in cases.items():
        a=c["data"]
        keep=("t","x","v","V","Q","p","pe","Fn","Fm","D","Dcage","Lring","wall","radial_gap","energy")
        matcases[name]={k:a[k] for k in keep}
    savemat(OUT/"ALBA_Step3_reference.mat",dict(parameters=p,cases=matcases),do_compression=True)
    from plot_step3 import make_plot
    make_plot(p,cases,fit)
    print(json.dumps({"resistance_fit":fit,"verification":verification},ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
