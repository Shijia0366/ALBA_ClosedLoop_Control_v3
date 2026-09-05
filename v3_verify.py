"""Independent frozen-Python reference execution. Never used as Simulink output."""
from pathlib import Path
import sys, json, math, time, hashlib, itertools
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'.python_deps'))
sys.path.insert(0,str(ROOT/'frozen_reference'))
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.io import loadmat, savemat
import stop_study as ss
PLANT=ss.Plant()

def write(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf-8')

class Command:
    def __init__(self,count,n=0):self.count,self.n=float(count),float(n)
    def scalar(self,t):return self.count,self.n

FIELDS=['phi','torque','Fm','Fn','fr','acc','tension','loss','slope','Qdot','pe','pg','A','D','y','dpdV','overlap','pressure']
UNITS=['rad','N*m','N','N','N','m/s^2','N','Pa','Pa*s/m^3','m^3/s^2','Pa','Pa','m^2','m','m','Pa/m^3','m','Pa']
ATOL=np.array([1e-10,1e-11,1e-8,1e-8,1e-10,1e-8,1e-8,1e-7,1e-2,1e-13,1e-7,1e-8,1e-12,1e-11,1e-11,1e-2,1e-11,1e-7])

def points():
    old=json.loads((ROOT/'stage9b_port/reference_points.json').read_text())
    inp=[[p[k] for k in ['x','v','V','Q','count','sealed']] for p in old['points']]
    # Signed overlap, just touching, actual force peak, low volume and phase boundary.
    for V,Q,v,gap,phi,sealed in itertools.product(
        [500e-6,300e-6,102e-6,PLANT.p['Vmin']+1e-10],
        [0.,3e-6,15e-6],[0.,.00311,PLANT.cfg.radius_m*2*math.pi-1e-9],
        [-.001,-1e-9,0.,1e-9,116.117/2e5,120/2e5],
        [0.,math.radians(78.6),math.pi/2-1e-8], [0,1]):
        yy=PLANT.geometry(V)[4];x=yy+gap
        count=(phi/50+x/.0025)*3200/(2*math.pi)
        inp.append([x,v,V,Q,count,sealed])
    # Source-order guard surfaces, including points just outside (point tests only).
    V=102e-6;yy=PLANT.geometry(V)[4]
    for j in range(7):
        for eps in [-1e-8,0,1e-8]:
            x=yy+116/2e5;v=.002;Q=3e-6;vv=V;phi=1.
            if j==0:phi=math.pi/2+eps
            elif j==1:x=PLANT.p['xmax']+eps
            elif j==2:vv=PLANT.p['Vmin']+eps*1e-3
            elif j==3:x=-1e-6+eps
            elif j==4:v=.0025*(2*math.pi+eps)
            elif j==5:phi=0.;x=yy-1e-4;v=eps
            else:Q=-1e-13+eps*1e-6
            count=(phi/50+x/.0025)*3200/(2*math.pi)
            inp.append([x,v,vv,Q,count,0])
    obs=[];rhs=[];jac=[];guards=[]
    for x,v,V,Q,count,sealed in inp:
        z=np.array([x,v,V,Q,0.,0.]);cmd=Command(count)
        o=PLANT.obs(0,z,cmd,sealed=bool(sealed))
        obs.append([o[k] for k in FIELDS]);rhs.append(PLANT.rhs(0,z,cmd,sealed=bool(sealed)))
        jac.append(PLANT.jac(0,z,cmd,sealed=bool(sealed)));guards.append(PLANT.guard_values(0,z,cmd,sealed=bool(sealed)))
    savemat(ROOT/'validation/points_v3.mat',dict(inputs=inp,obs=obs,rhs=rhs,jac=jac,guards=guards,obs_atol=ATOL))
    write(ROOT/'validation/point_test_spec.json',dict(points=len(inp),obs_fields=FIELDS,units=UNITS,atol=ATOL.tolist(),rtol=2e-10,
        scope='Includes inadmissible boundary-neighbour test points for port equivalence ONLY; not simulated beyond guards'))
    print('POINTS',len(inp),flush=True)

def reference(t,z,S,closed):
    r=S['target_mL']-S['tail_comp_mL']-z[4]*1e6
    qr=.5*S['q_plateau']*(1-math.cos(math.pi*np.clip(t/S['rise_s'],0,1)))
    u=np.clip((r-S['r_hold_mL'])/(S['R0_mL']-S['r_hold_mL']),0,1)
    qd=S['q_close']+.5*(S['q_plateau']-S['q_close'])*(1-math.cos(math.pi*u))
    return 0. if closed else min(qr,qd)

def control(t,z,S,closed):
    qref=reference(t,z,S,closed);Q=z[3]*1e6;n=z[6];ni=z[7]
    Vest=PLANT.p['V0']-z[4];A=PLANT.geometry(Vest)[1]
    nff=qref*1e-6/(2*math.pi*.0025/3200*A);e=qref-Q
    mode=int(S['mode']);raw=S['n_fixed'] if mode==1 else S['Kp']*e+ni+(nff if mode==3 else 0)
    term=False
    if S.get('terminal_manager',0):
        r=S['target_mL']-S['tail_comp_mL']-z[4]*1e6
        if mode==1:
            u=np.clip((r-S['r_hold_mL'])/(S['R0_mL']-S['r_hold_mL']),0,1)
            raw=S['n_fixed']*(S['q_close']+.5*(S['q_plateau']-S['q_close'])*(1-math.cos(math.pi*u)))/S['q_plateau']
        term=r<=S['stop_remaining_mL']
    des=0. if closed or term else raw;a=S['a_brake'] if closed else S['a_max']
    flag=float(des>S['n_max'])-float(des<0);des=np.clip(des,0,S['n_max'])
    want=(des-n)/S['tau_slew'];dn=float(np.clip(want,-a,a))
    if (n<=0 and dn<0) or (n>=S['n_max'] and dn>0):dn=0.
    dni=-S['Kaw']*ni if mode==1 or closed or term else S['Ki']*e+S['Kaw']*(n-raw)
    return dn,dni,qref,raw,e,flag,float(abs(want)>a+1e-9),float(term)

def tail(t,z,S,tc):
    if z[8]<=0 or z[9]>=S['tail_vol_mL']:return 0.
    T=2*S['tail_vol_mL']/z[8];u=np.clip(t-tc,0,T)
    return .5*z[8]*(1+math.cos(math.pi*u/T))

def rhs(t,z,S,closed,tc):
    phy=np.r_[z[:4],0.,0.];f=PLANT.rhs(t,phy,Command(z[5],z[6]),sealed=closed)
    dn,dni,*_=control(t,z,S,closed)
    return np.r_[f[:5],z[6],dn,dni,0. if closed else (z[3]*1e6-z[8])/S['tau_q0'],tail(t,z,S,tc) if closed else 0.]

def jacobian(t,z,S,closed,tc):
    J=np.zeros((10,10));phy=np.r_[z[:4],0.,0.];cmd=Command(z[5],z[6])
    j=PLANT.jac(t,phy,cmd,sealed=closed);J[:5,:4]=j[:5,:4]
    o=PLANT.obs(t,phy,cmd,sealed=closed)
    J[1,5]=PLANT.cfg.torque_cap*math.cos(o['phi'])/.0025*50*2*math.pi/3200/PLANT.M
    J[5,6]=1.
    for k,h in [(2,1e-12),(3,1e-12),(4,1e-12),(6,1e-3),(7,1e-3)]:
        zp=z.copy();zm=z.copy();zp[k]+=h;zm[k]-=h
        J[6:8,k]=(np.array(control(t,zp,S,closed)[:2])-np.array(control(t,zm,S,closed)[:2]))/(2*h)
    if not closed:J[8,3]=1e6/S['tau_q0'];J[8,8]=-1/S['tau_q0']
    return J

def guards(t,z,closed):
    return PLANT.guard_values(t,np.r_[z[:4],0.,0.],Command(z[5]),sealed=closed)[:6 if closed else 7]

def run_python(S,tmax=120,tight=False):
    rtol=2e-10 if tight else 2e-8;fac=.01 if tight else 1
    atol=np.array([2e-11,2e-10,2e-13,2e-14,2e-13,1e-6,1e-4,1e-4,1e-7,1e-9])*fac
    step=.005 if tight else .025;pieces=[];status='timeout';tc=math.nan
    z=np.array([0,0,500e-6,0,0,0,0,0,0,0.],dtype=float)
    def integrate(a,b,z,closed,closing=False):
        ev=[]
        for i in range(6 if closed else 7):
            def event(t,z,i=i):return guards(t,z,closed)[i]
            event.terminal=True;event.direction=-1;ev.append(event)
        if closing:
            def vol(t,z):return S['target_mL']-S['tail_comp_mL']-z[4]*1e6
            vol.terminal=True;vol.direction=-1;ev.append(vol)
        sol=solve_ivp(lambda t,z:rhs(t,z,S,closed,tc),(a,b),z,method='Radau',jac=lambda t,z:jacobian(t,z,S,closed,tc),
                      rtol=rtol,atol=atol,max_step=step,dense_output=True,events=ev)
        if not sol.success:raise RuntimeError(sol.message)
        pieces.append((closed,sol));return sol
    sol=integrate(0,tmax,z,False,True)
    if not len(sol.t_events[-1]):
        status='boundary_triggered' if any(len(v) for v in sol.t_events[:-1]) else 'timeout'
        return pieces,dict(status=status,t_close=tc,t_done=sol.t[-1])
    tthreshold=sol.t[-1]
    # Exact monotone v2 latch: comparator stays high, threshold .5 / gain.
    lag=.5/S['latch_gain'];sol=integrate(tthreshold,tthreshold+lag,sol.y[:,-1],False)
    if sol.t[-1]<tthreshold+lag-1e-10:return pieces,dict(status='boundary_triggered',t_close=tc,t_done=sol.t[-1])
    tc=sol.t[-1];z=sol.y[:,-1].copy();z[3]=0. # only Q resets, no mechanical reset
    T=2*S['tail_vol_mL']/z[8];tend=min(tmax,tc+max(T,2.5)+.5)
    sol=integrate(tc,tend,z,True)
    if any(len(v) for v in sol.t_events):return pieces,dict(status='boundary_triggered',t_close=tc,t_done=sol.t[-1])
    def quiet_margin(t):
        zz=sol.sol(t);o=PLANT.obs(t,np.r_[zz[:4],0.,0.],Command(zz[5]),sealed=True)
        return max(abs(zz[1]/.0025)/S['stop_omega_rad_s'],abs(o['acc'])/S['stop_acc_m_s2'])-1
    grid=np.unique(np.r_[sol.t,np.arange(tc,sol.t[-1],.001),sol.t[-1]])
    vals=np.array([quiet_margin(t) for t in grid]);cross=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
        if fa*fb<0:cross.append(brentq(quiet_margin,a,b,xtol=1e-11))
    bounds=[tc]+cross+[sol.t[-1]];tm=math.nan;done=math.nan
    tcmd=tc
    if sol.sol(tc)[6]>S['stop_n_steps_s']:
        tcmd=brentq(lambda t:sol.sol(t)[6]-S['stop_n_steps_s'],tc,sol.t[-1],xtol=1e-11)
    for a,b in zip(bounds[:-1],bounds[1:]):
        if quiet_margin((a+b)/2)<=0 and b-a>=S['stop_dwell_s']:
            candidate=max(max(a,tcmd)+S['stop_dwell_s'],tc+T)
            if candidate<=b:tm=a;done=candidate;break
    status='normal_completion' if np.isfinite(done) else 'timeout'
    if status=='normal_completion':
        # Distinct event times over the WHOLE trajectory, including a managed
        # stop before closure. No reset of mechanical states at the valve event.
        def state_at(t):
            for cc,pp in reversed(pieces):
                if t>=pp.t[0]-1e-12:return pp.sol(t),cc
            return pieces[0][1].sol(t),False
        def margin(t,mechanical):
            zz,cc=state_at(t)
            if not mechanical:return abs(zz[6])/S['stop_n_steps_s']-1
            oo=PLANT.obs(t,np.r_[zz[:4],0.,0.],Command(zz[5]),sealed=cc)
            return max(abs(zz[1]/.0025)/S['stop_omega_rad_s'],abs(oo['acc'])/S['stop_acc_m_s2'])-1
        tt=np.unique(np.r_[np.arange(0,done,.001),*[pp.t[pp.t<=done] for _,pp in pieces],done])
        event_times=[]
        for mechanical in [False,True]:
            mm=np.array([margin(t,mechanical) for t in tt]);bad=np.flatnonzero(mm>0)
            if len(bad) and bad[-1]<len(tt)-1:
                j=bad[-1];event_times.append(brentq(lambda t:margin(t,mechanical),tt[j],tt[j+1],xtol=1e-11))
            else:event_times.append(math.nan)
        tcmd,tm=event_times
    return pieces,dict(status=status,t_threshold=tthreshold,t_close=tc,t_tail_end=tc+T,t_command_stop=tcmd,t_mechanical_stop=tm,t_done=done if np.isfinite(done) else sol.t[-1])

def eval_pieces(pieces,tt):
    zz=np.empty((len(tt),10));closed=np.zeros(len(tt),bool)
    for c,sol in pieces:
        mask=(tt>=sol.t[0]-1e-12)&(tt<=sol.t[-1]+1e-12)
        zz[mask]=sol.sol(tt[mask]).T;closed[mask]=c
    return zz,closed

def channels(tt,zz,closed,S,tc):
    yy=np.zeros((len(tt),46));P=PLANT
    for k,(t,z,c) in enumerate(zip(tt,zz,closed)):
        o=P.obs(t,np.r_[z[:4],0.,0.],Command(z[5],z[6]),sealed=bool(c))
        dn,dni,qr,raw,e,sat,slew,term=control(t,z,S,c);qt=tail(t,z,S,tc) if c else 0
        g=np.array(guards(t,z,c));g=np.r_[g,np.inf] if c else g
        nv=sum(g<=0);gf=int(np.flatnonzero(g<=0)[0]+1) if nv else 0
        yy[k,:32]=[qr,z[3]*1e6,z[2]*1e6,500-z[2]*1e6,o['Fn'],z[6],raw,z[0],z[1],float(c),qt,z[9],500-z[2]*1e6+z[9],e,dn,z[7],z[8],tc if c else t,qt if c else z[3]*1e6,1. if c else 0.,float(500-z[2]*1e6>=S['target_mL']-S['tail_comp_mL']),sat,slew,o['phi'],o['torque'],z[1]/.0025,o['tension'],nv,gf,z[4]*1e6,z[5],o['acc']]
        yy[k,32:39]=g;yy[k,39]=z[0]-o['y'];yy[k,40]=o['pressure'];yy[k,41]=500-z[4]*1e6;yy[k,44]=term
    return yy

def trajectory(matfile):
    matfile=Path(matfile);d=loadmat(matfile,variable_names=['t','y','scenario','solver','diag'],simplify_cells=True)
    S=d['scenario'];tick=time.time();pieces,events=run_python(S,tmax=float(d.get('solver',{}).get('horizon_s',120)),tight=False)
    print('PYTHON',matfile.name,events,'wall',time.time()-tick,flush=True)
    end=events['t_done'];tt=np.unique(np.concatenate([sol.t for _,sol in pieces]));tt=tt[tt<=end];tt=np.r_[tt,end]
    zz,c=eval_pieces(pieces,tt);yy=channels(tt,zz,c,S,events['t_close'])
    outdir=ROOT/'validation'/'python_trajectories'/matfile.parent.name;outdir.mkdir(parents=True,exist_ok=True)
    savemat(outdir/matfile.name,dict(t=tt,y=yy,state=zz,scenario=S,events=events),do_compression=True)
    # Evaluate the independent dense solution at every Simulink time, side aware.
    ts=np.asarray(d['t']).ravel();ys=d['y'];valid=ts<=min(end,ts[-1]);ts=ts[valid];ys=ys[valid]
    zs,cs=eval_pieces(pieces,ts);yp=channels(ts,zs,cs,S,events['t_close'])
    same=(cs==(ys[:,9]>.5));delta=np.full_like(yp,np.nan)
    np.subtract(yp,ys,out=delta,where=np.isfinite(yp)&np.isfinite(ys))
    cols=[1,2,4,5,7,8,12,23,24,25,26,30,31,39,40]
    names=['Q_mL_s','V_mL','F_N','command_steps_s','x_m','v_m_s','collected_mL','phi_rad','torque_Nm','omega_rad_s','tension_N','count_STEP','acc_m_s2','signed_overlap_m','pressure_Pa']
    # Declared dimensional mixed absolute/relative test, not a unitless 1+|x| shortcut.
    atol=np.array([2e-3,2e-3,.03,.1,2e-7,2e-6,2e-3,2e-4,2e-5,8e-4,.03,.03,.01,2e-7,10.])
    phases={'startup':ts<=S['rise_s'],'plateau':(ts>S['rise_s'])&(ys[:,0]>=15-1e-6)&(~cs),
            'deceleration':(ts>S['rise_s'])&(ys[:,0]<15-1e-6)&(~cs),'postclose':cs}
    result={}
    for label,mask in phases.items():
        mask=mask&same
        if not np.any(mask):continue
        err=np.abs(delta[mask][:,cols]);lim=atol+2e-5*np.abs(yp[mask][:,cols])
        result[label]=dict(samples=int(sum(mask)),maximum_absolute=dict(zip(names,np.max(err,axis=0).tolist())),
            maximum_scaled_error=float(np.max(err/lim)),passed=bool(np.all(err<=lim)))
    sim_close=float(ts[np.flatnonzero(ys[:,9]>.5)[0]]) if np.any(ys[:,9]>.5) else math.nan
    rep=dict(events=events,simulink_close_s=sim_close,closure_time_difference_s=events['t_close']-sim_close,
             event_tolerance_s=.002,phase_results=result,absolute_tolerances=dict(zip(names,atol.tolist())),relative_tolerance=2e-5,
             different_branch_samples=int(sum(~same)),no_interpolation_across_Q_reset=True,
             passed=all(x['passed'] for x in result.values()) and
             ((not np.isfinite(events['t_close']) and not np.isfinite(sim_close)) or abs(events['t_close']-sim_close)<.002))
    # Event times are distinct. Mechanical event uses actual omega/acc only.
    diag=d['diag'];etable=[]
    for pk,sk in [('t_close','t_close_s'),('t_command_stop','command_stop_s'),('t_mechanical_stop','mechanical_stop_s'),('t_tail_end','tail_end_s'),('t_done','stop_time_s')]:
        if pk not in events:continue
        pv=float(events[pk]);sv=float(diag.get(sk,math.nan));bound=.002
        if pk=='t_mechanical_stop' and 'mechanical_event_bracket_s' in diag:
            bracket=np.asarray(diag['mechanical_event_bracket_s']);bound=max(bound,float(np.ptp(bracket))+.002)
        absent=not np.isfinite(pv) and not np.isfinite(sv)
        etable.append(dict(event=pk,python_s=pv,simulink_s=sv,difference_s=pv-sv,tolerance_s=bound,both_absent=absent,passed=bool(absent or abs(pv-sv)<=bound)))
    rep['event_comparison']=etable
    rep['simulink_status']=diag['status'];rep['same_status']=events['status']==diag['status']
    rep['passed']=rep['passed'] and rep['same_status'] and all(r['passed'] for r in etable)
    if events['status']=='boundary_triggered':
        zlast=pieces[-1][1].y[:,-1];gg=guards(end,zlast,pieces[-1][0])
        # Identify the zero-crossing event, NOT a minimum across different units.
        ev=pieces[-1][1].t_events
        triggered=[i+1 for i,v in enumerate(ev[:len(gg)]) if len(v)]
        rep['python_boundary']=dict(guard_indices=triggered,Q_actual_mL_s=zlast[3]*1e6,
            guard_margins=gg,guard_units=['rad','m','m^3','m','rad/s','N','m^3/s'][:len(gg)])
    write(outdir/(matfile.stem+'_comparison.json'),rep)
    savemat(outdir/(matfile.stem+'_aligned_errors.mat'),dict(t=ts,error=delta,python_y=yp,simulink_y=ys,same_branch=same),do_compression=True)
    print('TRAJECTORY PASS',rep['passed'],[(k,v['maximum_scaled_error']) for k,v in result.items()],flush=True)
    return rep['passed']

def audit():
    v2=ROOT/'provenance_v2'
    a=loadmat(ROOT/'ALBA_Stage8_Results.mat',simplify_cells=True);b=loadmat(v2/'ALBA_Stage8_Results.mat',simplify_cells=True)
    rows=[]
    for group in ['plant_parameters','hydraulic_parameters']:
        for name,value in a[group].items():
            if isinstance(value,(float,int,np.ndarray,np.number)):
                old=b[group][name];same=bool(np.array_equal(value,old));rows.append(dict(group=group,name=name,unchanged=same,new=np.asarray(value).tolist(),old=np.asarray(old).tolist()))
    for group in ['stiffness_coefficients','stiffness_powers','motor_radius_m','steps_per_rev']:
        rows.append(dict(group='ROM',name=group,unchanged=bool(np.array_equal(a[group],b[group])),new=np.asarray(a[group]).tolist(),old=np.asarray(b[group]).tolist()))
    hashes={}
    for f in ['alba9_const.m','stage9b_port/alba9b_const.m','ALBA_Stage8_Results.mat']:
        h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        hashes[f]=dict(v2=h(v2/f),v3=h(ROOT/f),unchanged=h(v2/f)==h(ROOT/f))
    assert all(r['unchanged'] for r in rows) and all(v['unchanged'] for v in hashes.values())
    write(ROOT/'validation/frozen_parameter_diff.json',dict(all_unchanged=True,rows=rows,hashes=hashes))
    import csv
    with (ROOT/'validation/frozen_parameter_diff.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=['group','name','unchanged','old','new']);w.writeheader();w.writerows(rows)
    print('FROZEN AUDIT PASS',len(rows),flush=True)

if __name__=='__main__':
    class Tee:
        def __init__(self,original,file):self.original,self.file=original,file
        def write(self,s):self.original.write(s);self.file.write(s);self.file.flush()
        def flush(self):self.original.flush();self.file.flush()
    logpath=ROOT/'logs'/('python_'+sys.argv[1]+'_'+time.strftime('%Y%m%dT%H%M%S')+'.log')
    logpath.parent.mkdir(exist_ok=True);logfile=logpath.open('w',encoding='utf-8')
    sys.stdout=Tee(sys.stdout,logfile);sys.stderr=Tee(sys.stderr,logfile)
    import scipy
    print('Python',sys.version,'NumPy',np.__version__,'SciPy',scipy.__version__,flush=True)
    cmd=sys.argv[1]
    if cmd=='points':points()
    elif cmd=='audit':audit()
    elif cmd=='trajectory':
        outcomes=[trajectory(arg) for arg in sys.argv[2:]]
        sys.exit(0 if all(outcomes) else 1)
