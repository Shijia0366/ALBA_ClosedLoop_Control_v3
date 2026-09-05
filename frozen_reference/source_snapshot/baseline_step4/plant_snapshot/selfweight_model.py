"""Conditional self-weight revision, SI units; NOT a validated folding model.

Measured fact: the entire ~0.270 kg cage is supported by the bladder and its
top button remains in contact. This does not determine a uniform pressure.
The diagnostic approximation is z_cg(V)-z_cg(V0)=beta*(Ds(V)-D0).
beta=1: concentrated dead load follows the top of a bottom-supported sphere.
beta=0.5: illustrative centroid-height sensitivity, NOT half the supported mass.
Actual z_cg(x,V,shape) is unmeasured; these are not rigorous bounds.

E_g=mg*beta*(Ds-D0), p_g=dE_g/dV=mg*beta/(pi*A), A=Ds**2/2.
Do not add mg to circumferential cable force. Gravity is counted exactly once.
Silicone pressure and geometry retain the OLD spherical assumption explicitly,
so post-collapse predictions expose model discrepancy, not a physical limit.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
import baseline


def parameters():
    p = baseline.parameters()
    p.update(cage_mass=0.270, support_fraction=1.0, gravity_height_slope=1.0,
             selfweight_model='conditional centroid-height law',
             initial_contact_fact='top button always contacts bladder',
             water_column_primed_assumption=True)
    return p


def geometry(V, p):
    a = baseline.shell(V, p)
    coefficient = p['cage_mass']*p['gravity']*p['gravity_height_slope']
    a['pg'] = coefficient/(np.pi*a['A'])
    a['Ug'] = coefficient*(a['D']-p['D0'])
    return a


def pressure_derivative(V, p):
    a = geometry(V, p)
    ri=(3*V/(4*np.pi))**(1/3)
    ro=a['D']/2
    ua=p['Dref_inner']/(2*ri); ub=p['Dref_outer']/(2*ro)
    dpe=(2*p['mu']/3)*(ua*(1+ua**3)/V-ub*(1+ub**3)/(V+p['Vsilicone']))
    dpg=-2*a['pg']/(3*(V+p['Vsilicone']))
    return dpe+dpg


def outputs(t, z, force, p):
    x,v,V=np.asarray(z)[:3]
    a=geometry(V,p)
    overlap=np.maximum(x-a['y'],0.)
    fn=p['contact_stiffness']*overlap
    pressure=a['pe']+a['pg']+fn/a['A']
    drive=pressure+p['rho']*p['gravity']*p['head']
    Q=np.maximum(drive,0.)/p['Rh']
    fm=baseline.force_input(t,force,p)
    fr=p['friction']*np.tanh(v/p['friction_velocity'])
    a.update(x=x,v=v,V=V,Fn=fn,p=pressure,drive=drive,Q=Q,Fm=fm,
             Ffriction=fr,overlap=overlap,radial_gap=(a['y']-x)/(2*np.pi),
             Dcage=p['D0']-x/np.pi,Lring=p['L0']-x)
    return a


def rhs(t,z,force,p):
    a=outputs(t,z,force,p)
    acceleration=(a['Fm']-p['damping']*a['v']-a['Ffriction']-a['Fn'])/p['mass_equivalent']
    dz=[a['v'],acceleration,-a['Q']]
    if len(z)==5:
        power=(a['Fm']*a['v']+p['rho']*p['gravity']*p['head']*a['Q']
               -p['damping']*a['v']**2-a['Ffriction']*a['v']-p['Rh']*a['Q']**2)
        dz += [a['Q'],power]
    return dz


def jacobian(t,z,force,p):
    a=outputs(t,z,force,p)
    V=float(z[2]); v=float(z[1]); area=float(a['A'])
    k=p['contact_stiffness'] if a['overlap']>0 else 0.
    fnV=k/area; darea=2/(np.pi*a['D'])
    qx=qV=0.
    if a['drive']>0:
        qx=k/(area*p['Rh'])
        qV=(pressure_derivative(V,p)+fnV/area-a['Fn']*darea/area**2)/p['Rh']
    df=p['friction']/p['friction_velocity']*(1-np.tanh(v/p['friction_velocity'])**2)
    J=np.zeros((len(z),len(z)))
    J[0,1]=1
    J[1,:3]=[-k/p['mass_equivalent'],-(p['damping']+df)/p['mass_equivalent'],-fnV/p['mass_equivalent']]
    J[2,:3]=[-qx,0,-qV]
    if len(z)==5:
        J[3,:3]=[qx,0,qV]
        h=p['rho']*p['gravity']*p['head']-2*p['Rh']*a['Q']
        J[4,:3]=[h*qx,a['Fm']-a['Ffriction']-v*df-2*p['damping']*v,h*qV]
    return J


def run_case(force,p,horizon):
    def limit(t,z): return p['xmax']-z[0]
    def sphere_boundary(t,z): return z[2]-p['Vmin']
    def extension(t,z): return z[0]+1e-6
    for e in (limit,sphere_boundary,extension):
        e.terminal=True; e.direction=-1
    sol=solve_ivp(lambda t,z:rhs(t,z,force,p),(0,horizon),[0.,0.,p['V0'],0.,0.],
                  method='Radau',rtol=p['rtol'],atol=p['atol'],dense_output=True,
                  max_step=p['max_step'],jac=lambda t,z:jacobian(t,z,force,p),
                  events=(limit,sphere_boundary,extension))
    if not sol.success: raise RuntimeError(sol.message)
    t=np.unique(np.r_[np.linspace(0,min(3,sol.t[-1]),501),sol.t,
                       np.arange(3,sol.t[-1],.5),sol.t[-1]])
    z=sol.sol(t); a=outputs(t,z,force,p)
    energy=(.5*p['mass_equivalent']*z[1]**2+baseline.shell_energy(z[2],p)
            +a['Ug']+.5*p['contact_stiffness']*a['overlap']**2)
    a.update(t=t,z=z,energy=energy)
    return dict(sol=sol,data=a,checks=dict(
        mass_balance_error_mL=float(np.max(np.abs(z[3]-(p['V0']-z[2])))*1e6),
        energy_balance_error_J=float(np.max(np.abs(energy-energy[0]-z[4]))),
        finite=bool(all(np.all(np.isfinite(a[k])) for k in ('x','v','V','p','Q'))),
        ended_at_s=float(sol.t[-1]),below_natural=bool(np.any(z[2]<p['Vnatural'])),
        physical_postcollapse_prediction_valid=False))


def evaluate(result,trial,p):
    sol=result['sol']; t=np.asarray(trial['elapsed_s_points'])
    if sol.t[-1]<max(t)-1e-6: raise RuntimeError('Model stopped before trial endpoint')
    V=sol.sol(t)[2]; pred=(p['V0']-V)*1e6; obs=np.asarray(trial['volume_mL'])
    err=pred-obs
    return dict(id=trial['id'],mass_g=trial['mass_g'],repeat=trial['repeat'],
                elapsed_s=t.tolist(),observed_mL=obs.tolist(),predicted_mL=pred.tolist(),
                final_observed_mL=float(obs[-1]),final_predicted_mL=float(pred[-1]),
                rmse_mL=float(np.sqrt(np.mean(err[1:]**2))),
                endpoint_error_mL=float(err[-1]))
