import re, sys, time
tex = open('/home/mbustamante/Research/magnus/resources/paper/main.tex').read()
code = re.search(r'\\begin\{lstlisting\}\[[^\]]*label=\{lst:validation\}[^\]]*\]\n(.*?)\\end\{lstlisting\}', tex, re.S).group(1)
for label, c in (('1e-12', code), ('1e-3 ', code.replace('rtol=1e-12, atol=1e-14,', 'rtol=1e-3, atol=1e-3,'))):
    parts = re.split(r'\n(?=P_\dnu = )', c)
    g = {}; exec(parts[0], g)
    for p in parts[1:]:
        exec(p.replace("energies(", "energies2(").replace('P_', 'W_', 1), dict(g, energies2=lambda lo, hi: g['energies'](lo, hi)[:2]))
    T = {}
    for p in parts[1:]:
        t = time.perf_counter(); exec(p, g); T[p.split(' =')[0]] = time.perf_counter() - t
    print(label, ' '.join('%s=%.3f' % kv for kv in T.items()), 'total=%.3f' % sum(T.values()))
