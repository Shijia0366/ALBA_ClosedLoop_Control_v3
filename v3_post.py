"""Side-aware numerical comparisons and unclipped scientific figures."""
from pathlib import Path
import sys,json,math
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'.python_deps'))
import numpy as np
from scipy.io import loadmat,savemat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})

def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def load(p):
    d=loadmat(p,variable_names=['t','y','scenario','diag'],simplify_cells=True)
    d['t']=np.asarray(d['t']).ravel()
    if d['y'].shape[0]!=len(d['t']):d['y']=d['y'].T
    return d
def interp_side(t,y,tt,closed):
    use=y[:,9]>.5 if closed else y[:,9]<.5
    t,y=t[use],y[use];unique,idx=np.unique(t[::-1],return_index=True);y=y[::-1][idx]
    return np.column_stack([np.interp(tt,unique,y[:,j]) for j in range(y.shape[1])])

def compare(a,b,label,kind):
    A=load(a);B=load(b);ta,ya=A['t'],A['y'];tb,yb=B['t'],B['y'];rise=float(A['scenario']['rise_s'])
    ca=ta[ya[:,9]>.5][0] if np.any(ya[:,9]>.5) else ta[-1]
    cb=tb[yb[:,9]>.5][0] if np.any(yb[:,9]>.5) else tb[-1]
    out=ROOT/'validation'/label;out.mkdir(parents=True,exist_ok=True)
    ts=[];aa=[];bb=[];ph=[]
    for closed in [False,True]:
        lo=max(ca,cb) if closed else 0.;hi=min(ta[-1],tb[-1]) if closed else min(ca,cb)-1e-9
        if hi<=lo:continue
        tt=np.unique(np.r_[np.arange(lo,hi,.0002),ta[(ta>=lo)&(ta<=hi)],tb[(tb>=lo)&(tb<=hi)],hi])
        x=interp_side(ta,ya,tt,closed);z=interp_side(tb,yb,tt,closed)
        phase=np.where(tt<=rise,0,np.where(x[:,0]>=15-1e-6,1,2));phase[:]=3 if closed else phase
        ts.append(tt);aa.append(x);bb.append(z);ph.append(phase)
    t=np.concatenate(ts);x=np.concatenate(aa);z=np.concatenate(bb);phase=np.concatenate(ph)
    cols=[1,2,4,5] if kind=='rom' else [1,2,4,5,7,8,23,24,25,31,39,40]
    names=['Q_mL_s','V_mL','F_N','command_STEP_s'] if kind=='rom' else ['Q_mL_s','V_mL','F_N','command_STEP_s','x_m','v_m_s','phi_rad','torque_Nm','omega_rad_s','acc_m_s2','overlap_m','pressure_Pa']
    err=x[:,cols]-z[:,cols];results={}
    tol=np.array([.002,.002,.03,.1,2e-7,2e-6,2e-4,2e-5,8e-4,.01,2e-7,10.]) if kind!='rom' else None
    for p,nam in enumerate(['startup','plateau','deceleration','postclose']):
        m=phase==p
        if not np.any(m):continue
        results[nam]={'maximum_absolute':dict(zip(names,np.max(abs(err[m]),axis=0).tolist()))}
        if tol is not None:
            ratio=abs(err[m])/(tol+2e-5*abs(z[m][:,cols]));results[nam]['maximum_scaled_error']=float(np.max(ratio));results[nam]['passed']=bool(np.all(ratio<=1))
    rep={'kind':kind,'a':str(a),'b':str(b),'closure_time_difference_s':float(ca-cb),'phases':results,
         'excluded_different_mode_window_s':[float(min(ca,cb)),float(max(ca,cb))], 'interpolation':'linear within each mode only; no Q reset crossing'}
    if tol is not None:rep['passed']=all(v['passed'] for v in results.values()) and abs(ca-cb)<.002
    write(out/'comparison.json',rep);savemat(out/'aligned_errors.mat',dict(t=t,errors=err,phase=phase),do_compression=True)
    fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True)
    for j,ax in enumerate(axs):
        col=0 if j==0 else 2
        for p in range(4):
            m=phase==p;ax.plot(t[m],err[m,col],lw=.8,color='#245c8a')
        ax.axhline(0,color='gray',lw=.5);ax.grid(alpha=.2);ax.set_ylabel(['Aligned Q error (mL/s)','Aligned contact-force error (N)'][j]);ax.axvline(rise,color='gray',ls=':')
    axs[-1].set_xlabel('Time (s)');fig.suptitle('ROM minus full model (not hardware)' if kind=='rom' else 'Nominal minus tightened full-model simulation')
    fig.tight_layout();fig.savefig(out/'aligned_errors.png',dpi=180);axs[-1].set_xlim(0,3);fig.savefig(out/'startup_aligned_errors.png',dpi=180);plt.close(fig)
    print(label,rep.get('passed','ROM discrepancy quantified'),results['startup'],flush=True)
    return rep

def plot_case(path):
    d=load(path);t,y,S=d['t'],d['y'],d['scenario'];out=path.parent/(path.stem+'_figures');out.mkdir(exist_ok=True)
    close=t[y[:,9]>.5][0] if np.any(y[:,9]>.5) else None
    for scope in ['full','startup','closure','braking']:
        if scope in ['closure','braking'] and close is None:continue
        fig,axs=plt.subplots(4,2,figsize=(12,11),sharex=True);a=axs.ravel()
        # Keep both sides of valve transition in the raw log. No smoothing.
        a[0].plot(t,y[:,0],ls='--',label='Q_ref');a[0].plot(t,y[:,1],label='Q_actual through valve');a[0].plot(t,y[:,18],ls=':',label='Final outlet')
        a[0].set_ylabel('Flow (mL/s)');a[0].legend(fontsize=7)
        a[1].plot(t,y[:,3],label='Bladder discharge');a[1].plot(t,y[:,12],label='Final collected');a[1].axhline(S['target_mL'],ls='--',color='gray');a[1].set_ylabel('Volume (mL)');a[1].legend(fontsize=7)
        a[2].plot(t,y[:,5]*60/3200,label='Command equivalent rpm');a[2].plot(t,y[:,25]*60/(2*np.pi),label='Actual rotor rpm');a[2].set_ylabel('Speed (rpm)');a[2].legend(fontsize=7)
        a[3].plot(t,y[:,24],label='Actual torque');a[3].axhline(.29658517619749203,ls='--',color='gray',label='Conditional envelope');a[3].plot(t,.29658517619749203-abs(y[:,24]),label='Torque headroom');a[3].set_ylabel('Torque (N m)');a[3].legend(fontsize=7)
        a[4].plot(t,y[:,4],label='Contact generalized force');a[4].plot(t,y[:,26],ls=':',label='Cable tension');a[4].set_ylabel('Force (N)');a[4].legend(fontsize=7)
        a[5].plot(t,y[:,23]*180/np.pi);a[5].axhline(90,ls='--',color='gray');a[5].set_ylabel('Electrical phase (deg)')
        a[6].plot(t,y[:,39]*1e3);a[6].axhline(0,ls='--',color='gray');a[6].set_ylabel('Signed overlap x-y (mm)')
        a[7].plot(t,y[:,40]/1000);a[7].set_ylabel('Bladder gauge pressure (kPa)')
        for ax in a:
            ax.grid(alpha=.2)
            if close is not None:ax.axvline(close,color='#aa5544',lw=.7,ls=':')
            if scope=='startup':ax.set_xlim(0,3)
            elif scope=='closure':ax.set_xlim(max(0,close-2),min(t[-1],close+2.5))
            elif scope=='braking':
                lo=max(0,close-.10);hi=min(t[-1],close+.30);ax.set_xlim(lo,hi)
                vals=[]
                for line in ax.lines:
                    xx=np.asarray(line.get_xdata());vv=np.asarray(line.get_ydata())
                    if len(xx)==len(t):vals.extend(vv[(xx>=lo)&(xx<=hi)].tolist())
                if vals:
                    low,high=min(vals),max(vals);pad=max((high-low)*.15,abs(high)*.02,1e-7)
                    ax.set_ylim(low-pad,high+pad)
        a[-1].set_xlabel('Time (s)');a[-2].set_xlabel('Time (s)')
        status=d.get('diag',{}).get('status','unknown')
        fig.suptitle(f"{path.parent.name} | {int(S['target_mL'])} mL | mode {int(S['mode'])} | {status}\nConditional model only; bench effective resistance unchanged",fontsize=11)
        fig.tight_layout(rect=(0,0,1,.95));fig.savefig(out/(scope+'.png'),dpi=160);plt.close(fig)

def poster(tag,mode):
    fig,axs=plt.subplots(2,1,figsize=(10,6.5),sharex=True)
    colors=['#138b8d','#d68832','#6e62ad']
    for target,color in zip([200,300,400],colors):
        d=load(ROOT/'results'/tag/f'case_{target}mL_mode{mode}.mat');t,y=d['t'],d['y']
        axs[0].plot(t,y[:,18],color=color,lw=1.7,label=f'{target} mL target')
        axs[0].plot(t,y[:,0],color=color,lw=.7,ls='--',alpha=.5)
        axs[1].plot(t,y[:,12],color=color,lw=1.8);axs[1].axhline(target,color=color,lw=.7,ls=':')
    axs[0].set_ylabel('Final outlet flow (mL/s)');axs[0].legend(ncol=3,frameon=False)
    axs[1].set_ylabel('Collected volume (mL)');axs[1].set_xlabel('Time (s)')
    for ax in axs:ax.grid(alpha=.2)
    fig.suptitle('ALBA | nominal-model closed-loop simulation',fontsize=15)
    fig.text(.5,.071,'Solid: final-outlet flow  |  Dashed: Q_ref  |  Dotted: collection target',ha='center',fontsize=8)
    fig.text(.5,.025,'15 mL/s: Upper-feasible / near-limit case under nominal model assumptions.\nNot hardware validation; terminal volume follows ideal valve and assumed 2 mL tail.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.94));out=ROOT/'poster';out.mkdir(exist_ok=True)
    for ext in ['png','svg','pdf']:fig.savefig(out/f'ALBA_final_flow_volume.{ext}',dpi=300)
    plt.close(fig)

if __name__=='__main__':
    if sys.argv[1]=='compare':compare(Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4],sys.argv[5])
    elif sys.argv[1]=='plots':
        for arg in sys.argv[2:]:
            for path in (ROOT/'results'/arg).glob('case*.mat'):
                if 'failure' not in path.name:plot_case(path)
    elif sys.argv[1]=='poster':poster(sys.argv[2],int(sys.argv[3]))
    elif sys.argv[1]=='checks':
        for mode in [1,2,3]:
            file=f'case_300mL_mode{mode}.mat'
            tight=ROOT/'results/D_convergence_tight'/file
            if tight.exists():compare(ROOT/'results/D_diagnostic_baseline'/file,tight,f'convergence_300_mode{mode}','tight')
