"""Regenerate summary CSVs and standalone figures from the saved experiments."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .methods import generate, measure, size_counts, fit_counts

BASE=Path(__file__).resolve().parent
OUT=BASE/'results'


def wilson(hits,n):
    z=1.959964;p=hits/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half


def main():
    fits=pd.read_csv(OUT/'fits.csv')
    primary=fits[(fits.fit_min==1)&(fits.fit_max==30)].copy()
    primary['mle_and_gof']=primary.signature_mle & (primary.ks_iid_p>.05)
    primary['gof_reject']=primary.ks_iid_p<=.05
    summaries=[]
    for keys,g in primary.groupby(['hours','condition','policy']):
        r=dict(zip(['hours','condition','policy'],keys),n=len(g))
        for col in ['sigma_single','tau_regression','tau_mle','ks','aic_powerlaw_minus_lognormal']:
            r[col+'_mean']=g[col].mean()
            r[col+'_sd']=g[col].std(ddof=1)
            r[col+'_min']=g[col].min();r[col+'_max']=g[col].max()
        for col in ['signature_regression','signature_mle','mle_and_gof','gof_reject']:
            r[col+'_hits']=int(g[col].sum())
            r[col+'_fraction']=float(g[col].mean())
            r[col+'_ci_low'],r[col+'_ci_high']=wilson(g[col].sum(),len(g))
        summaries.append(r)
    pd.DataFrame(summaries).to_csv(OUT/'decision_summary.csv',index=False)
    intervals=pd.read_csv(OUT/'interval_tests.csv')
    rows=[]
    for (condition,width),g in intervals.groupby(['condition','window_ms']):
        hits=int((g.rank_p<=.05).sum());low,high=wilson(hits,len(g))
        rows.append(dict(condition=condition,window_ms=width,hits=hits,n=len(g),
                         fraction=hits/len(g),ci_low=low,ci_high=high))
    interval_summary=pd.DataFrame(rows)
    interval_summary.to_csv(OUT/'interval_summary.csv',index=False)

    cal=pd.read_csv(OUT/'calibration.csv')
    cal.assign(reject=cal.ks_iid_p<=.05).groupby(['distribution','n']).agg(
        runs=('seed','size'),rejects=('reject','sum'),tau_mean=('tau_mle','mean')
    ).to_csv(OUT/'calibration_summary.csv')

    # A compact, reusable research figure. Error bars are Wilson intervals
    # over independent runs, not uncertainty over individual avalanches.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12.5,5.2),layout='constrained')
    policies=['iei_rounded','fixed_10','selected_mle']
    labels=['Prescribed IEI','Fixed 10 ms','Bin selected\nusing MLE']
    colors=['#7756a6','#178386','#c47a28']
    for k,(col,label) in enumerate([('signature_regression','Regression + LLR'),
                                   ('signature_mle','MLE + LLR'),
                                   ('mle_and_gof','MLE + LLR + fit check')]):
        fracs=[];lowers=[];uppers=[]
        for policy in policies:
            g=primary[(primary.hours==3)&(primary.condition=='mixed')&(primary.policy==policy)]
            hits=g[col].sum();low,high=wilson(hits,len(g));p=hits/len(g)
            fracs.append(p);lowers.append(max(0,p-low));uppers.append(max(0,high-p))
        axes[0].bar(np.arange(3)+(k-1)*.24,fracs,.23,label=label,color=colors[k],
                    yerr=[lowers,uppers],capsize=3,error_kw={'linewidth':.8})
    axes[0].set(xticks=np.arange(3),xticklabels=labels,ylim=(0,1.1),ylabel='Fraction meeting descriptive screen',
                title='A  Mixed-drive null: 20 independent 3-hour runs')
    axes[0].legend(loc='upper center',bbox_to_anchor=(.5,-.14),frameon=False,fontsize=9)
    cs=['homogeneous','uniform','mixed','capped','refractory_1.0']
    lab=['Homogeneous\nnull','Uniform-drive\nnull','Mixed-drive\nnull','Capped\nbranching','Refractory\nbranching']
    for k,width in enumerate([20,40,100]):
        sub=interval_summary[interval_summary.window_ms==width].set_index('condition').loc[cs]
        axes[1].bar(np.arange(len(cs))+(k-1)*.24,sub.fraction,.23,label=f'{width} ms windows',color=colors[k],
                    yerr=np.maximum(0,np.array([sub.fraction-sub.ci_low,sub.ci_high-sub.fraction])),
                    capsize=2,error_kw={'linewidth':.8})
    axes[1].set(xticks=np.arange(len(cs)),xticklabels=lab,ylim=(0,1.1),ylabel='Fraction rejecting timing exchangeability',
                title='B  Interval jitter: 20 independent 15-minute runs')
    axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.14),frameon=False,fontsize=9)
    for ax in axes: ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    for ext in ['png','svg']:fig.savefig(OUT/f'decision_summary.{ext}',dpi=200)
    plt.close(fig)

    # Recreate the headline recording solely for distribution plotting.
    _,times,_=generate('mixed',70.,42)
    f=measure(times,10);counts,_=size_counts(f['size'].to_numpy(),1,30)
    fit,lp=fit_counts(counts);x=np.arange(1,31);observed=counts/counts.sum()
    long=pd.read_csv(OUT/'long_fits.csv')
    row=long[(long.seed==42)&(long.policy=='fixed_10')&(long.fit_min==1)&(long.fit_max==30)].iloc[0]
    lr=-row.lognormal_a*np.log(x)-row.lognormal_b*np.log(x)**2
    ln=np.exp(lr-lr.max());ln/=ln.sum()
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.4),layout='constrained')
    axes[0].loglog(x,observed,'o',ms=4,color='#303744',label='Observed mixed-drive null')
    axes[0].loglog(x,np.exp(lp),color='#178386',label=f'Power-law MLE: tau={fit["tau_mle"]:.3f}')
    axes[0].loglog(x,ln,'--',color='#c47a28',label='Lognormal-shaped fit')
    axes[0].set(xlabel='Avalanche size',ylabel='Probability conditional on sizes 1–30',title='A  70 hours, seed 42, 10 ms bins')
    axes[0].legend(frameon=False,fontsize=8)
    sub=long[(long.policy=='fixed_10')&(long.fit_min==1)&(long.fit_max==30)]
    for i,(col,label,color) in enumerate([('tau_regression','Log-log regression','#7756a6'),('tau_mle','Maximum likelihood','#178386')]):
        axes[1].scatter(i+np.linspace(-.08,.08,len(sub)),sub[col],color=color,s=32,label=label)
    axes[1].axhline(1.5,ls=':',color='#555555',label='Reference exponent 1.5')
    axes[1].set(xticks=[0,1],xticklabels=['Regression','MLE'],ylabel='Estimated exponent',
                title='B  Six independent 70-hour null recordings',xlim=(-.5,1.5),ylim=(1.2,1.6))
    axes[1].grid(axis='y',alpha=.2)
    for ext in ['png','svg']:fig.savefig(OUT/f'estimator_gap.{ext}',dpi=200)
    plt.close(fig)
    print('Summary CSVs and four figure files written.')


if __name__=='__main__':main()
