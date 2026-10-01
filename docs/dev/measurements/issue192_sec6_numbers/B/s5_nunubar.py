import numpy as np, warnings, time
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob
import magnus.earth as earth
import magnus.globaldefs as gd
costhz = -0.9
L = earth.distance_traveled_inside_earth(costhz) \
 *gd.UNIT_KM
kw = dict(costhz=costhz, L=L, nu_i=gd.NUMU,
          nu_f=gd.NUMU, rtol=1.0e-8, atol=1.0e-10)
def both_signs(wrapper, E, **params):
    return [wrapper(E, nubar=nubar, **params, **kw)
            for nubar in (False, True)]
osc = gd.load_nufit_params('NuFIT 6.1')
E = np.logspace(0, np.log10(40), 260)*gd.UNIT_GEV
t=time.time()
P2, P2bar = both_signs(oscprob.osc_prob_2nu_earth, E,
                       sth=osc['s13'], Dm2=osc['D31'])
P3, P3bar = both_signs(oscprob.osc_prob_3nu_earth, E)
E_TEV = np.logspace(0, np.log10(30), 200)*gd.UNIT_TEV
s4 = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0,
          D41=1.0)
s5 = dict(**s4, s15=np.sqrt(0.06), s25=np.sqrt(0.06),
          s35=0.0, D51=1.7)
P4, P4bar = both_signs(oscprob.osc_prob_4nu_earth,
                       E_TEV, **s4)
P5, P5bar = both_signs(oscprob.osc_prob_5nu_earth,
                       E_TEV, **s5)
print('ran in', time.time()-t)
np.savez('nunubar.npz', E=E/gd.UNIT_GEV, ET=E_TEV/gd.UNIT_TEV, P2=P2,P2bar=P2bar,P3=P3,P3bar=P3bar,P4=P4,P4bar=P4bar,P5=P5,P5bar=P5bar)
Eg = E/gd.UNIT_GEV; Et = E_TEV/gd.UNIT_TEV
P2,P2bar,P3,P3bar,P4,P4bar,P5,P5bar=map(np.asarray,(P2,P2bar,P3,P3bar,P4,P4bar,P5,P5bar))
k=np.argmin(P2); print('2nu numu min %.4f at %.3f GeV; nubar min %.4f'%(P2[k],Eg[k],P2bar.min()))
d=np.abs(P3-P3bar); m=Eg<5; k=np.argmax(d*m); print('3nu max diff below 5 GeV %.4f at %.3f'%(d[k],Eg[k]))
for th in (10,15,20,25,30): print(' max diff above',th,'GeV: %.4f'%d[Eg>th].max())
k=np.argmin(P4bar); print('3+1 nubar min %.4f at %.3f TeV; nu min %.4f'%(P4bar[k],Et[k],P4.min()))
# 3+2 dips: local minima of P5bar
idx=[i for i in range(1,len(Et)-1) if P5bar[i]<P5bar[i-1] and P5bar[i]<P5bar[i+1]]
print('3+2 nubar local minima', [(round(Et[i],3), round(P5bar[i],4)) for i in idx])
print('3+2 nu min %.4f'%P5.min())
for th in (5,10,15,20): print(' >%d TeV max|P4-P4bar| %.4f  max|P5-P5bar| %.4f'%(th, np.abs(P4-P4bar)[Et>th].max(), np.abs(P5-P5bar)[Et>th].max()))
# widths (FWHM-like: below midpoint between min and 1)
def width(E,P,lo,hi):
    m=(E>lo)&(E<hi); Pm=P[m]; Em=E[m]; k=np.argmin(Pm); half=(Pm[k]+1)/2
    i=k
    while i>0 and Pm[i]<half: i-=1
    j=k
    while j<len(Pm)-1 and Pm[j]<half: j+=1
    return Em[k],Pm[k],Em[i],Em[j],np.log10(Em[j]/Em[i])
print('3+1 dip width', width(Et,P4bar,0.5,30))
