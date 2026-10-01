import numpy as np, warnings, time
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.avgprob as avgprob
from shockdef import *
OSC = gd.load_nufit_params('NuFIT 6.1')
E_SN = np.logspace(np.log10(5.0), np.log10(60.0), 60)*gd.UNIT_MEV
KW = dict(L0=L0, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True)
none = np.asarray(oscprob.osc_prob_matter_std_potential(3, undisturbed_ne, E_SN, L1, OSC, average=True, **KW), float)
print('no shock P in %.4f..%.4f' % (none.min(), none.max()))
out={'none':none}
for key, width in (('0p07',1e-6),('70',1e-3)):
    ne=sn_shock_ne(width); bp=np.array([L0]+sorted(edges_for(width*7e4))+[L1])
    f=lambda e: np.asarray(oscprob.osc_prob_matter_std_potential(3, ne, e, L1, OSC, t_breakpoints=bp, **KW))
    t=time.time(); ms=[];ss=[]
    for e in E_SN:
        m,s=avgprob.averaged_probabilities_numerically(f, float(e)); ms.append(float(m)); ss.append(float(s))
    out[key]=np.array(ms); out[key+'_sem']=np.array(ss)
    print(key, '%.0fs' % (time.time()-t), 'mean %.3f..%.3f sem max %.3f' % (min(ms),max(ms),max(ss)))
np.savez('energy.npz', E=E_SN/gd.UNIT_MEV, **out)
Em=E_SN/gd.UNIT_MEV
for i in range(0,60,3):
    print('%6.2f  none %.3f  sharp %.3f  wide %.3f' % (Em[i], none[i], out['0p07'][i], out['70'][i]))
