import re, json, warnings, inspect, numpy as np
warnings.simplefilter('ignore')
import magnus.oscprob as op
print('osc_prob default magnus_exp_order =', inspect.signature(op.osc_prob).parameters['magnus_exp_order'].default)
cache = json.load(open('/home/mbustamante/Research/magnus/notebooks/paper_figure_cache.json'))['oracle']
val = cache['value'] if 'value' in cache else cache
print('cache entry measured:', cache.get('measured'), 'on', cache.get('machine'))
tex = open('/home/mbustamante/Research/magnus/resources/paper/main.tex').read()
code = re.search(r'\\begin\{lstlisting\}\[[^\]]*label=\{lst:validation\}[^\]]*\]\n(.*?)\\end\{lstlisting\}', tex, re.S).group(1)
g = {}; exec(code, g)
for d, name in ((2, 'P_2nu'), (3, 'P_3nu'), (4, 'P_4nu'), (5, 'P_5nu')):
    v = val[str(d)]
    Ef = np.array(v['E_plot']); Pf = np.array(v['P_plot'])[:, 0, 0]
    El = g['energies'](*{2: (0.0005, 0.05), 3: (0.002, 0.2), 4: (2.0, 20.0), 5: (2.0, 20.0)}[d])
    Pl = np.asarray(g[name]).ravel()
    assert np.allclose(Ef, El, rtol=1e-14), d
    print('d=%d  max|listing - Fig.1 curve| = %.1e   Fig.1 residual max %.1e   oracle floor max %.1e'
          % (d, np.max(np.abs(Pl - Pf)), max(v['resid']), max(v['floor'])))
