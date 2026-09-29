"""Run Listing 1 exactly as printed in main.tex, timing each curve; optional tolerance override."""
import re, sys, time
tex = open('/home/mbustamante/Research/magnus/resources/paper/main.tex').read()
m = re.search(r'\\begin\{lstlisting\}\[[^\]]*label=\{lst:validation\}[^\]]*\]\n(.*?)\\end\{lstlisting\}', tex, re.S)
code = m.group(1)
if len(sys.argv) > 1 and sys.argv[1] == 'default':
    code = code.replace('rtol=1e-12, atol=1e-14,', 'rtol=1e-3, atol=1e-3,')
    assert 'rtol=1e-3' in code
parts = re.split(r'\n(?=P_\dnu = )', code)
g, T = {}, {}
exec(parts[0], g)
for p in parts[1:]:
    name = p.split(' =')[0]
    t = time.perf_counter(); exec(p, g); T[name] = time.perf_counter() - t
print(' '.join('%s=%.3f' % kv for kv in T.items()), 'total=%.3f' % sum(T.values()))
