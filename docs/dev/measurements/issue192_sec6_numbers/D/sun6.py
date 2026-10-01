from sun3 import *
bs=np.arange(0.70,0.95,0.002); pv=np.array([Pvary(300.0,b) for b in bs])
i=pv.argmin(); print('averaged_varying 300 GeV: min %.4f at b=%.3f; at 0.81: %.4f' % (pv[i], bs[i], np.interp(0.81,bs,pv)))
