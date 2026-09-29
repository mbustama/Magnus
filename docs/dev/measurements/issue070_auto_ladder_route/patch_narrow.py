p = 'src/magnus/oscprob.py'; s = open(p).read()

def splice(start, end, new):
    global s
    i0 = s.index(start); i1 = s.index(end, i0)
    s = s[:i0] + new + s[i1:]

D = '/tmp/claude-1000/-home-mbustamante-Research-magnus/a0d2e3ee-05d1-4f5e-8936-3df18ac27060/scratchpad/i70/'
splice('AUTO_LADDER_MAX_PHASE = 1.0e4\n', 'AUTO_LADDER_MIN_TOLERANCE = 1.0e-6\n', open(D + 'blk_const.txt').read())
splice('def _estimated_phase(', 'def _osc_prob_hybrid_dispatch(', open(D + 'blk_estimate2.txt').read())

old = '''    prefer = (_auto_prefers_ladder(H_at_energy, energy_arr, L_arr, L0, rtol, atol)
              if strategy == 'auto' else None)
'''
assert s.count(old) == 2
i = s.index(old)
s = s[:i] + '''    prefer = (_auto_prefers_ladder(H_at_energy, energy_arr, L_arr, L0, rtol, atol,
                  _resolve_max_n_slabs(scan_kwargs.get('max_n_slabs'), integration_method))
              if strategy == 'auto' else None)
''' + s[i + len(old):]
i = s.index(old)
s = s[:i] + '''    prefer = (_auto_prefers_ladder(H_at_energy, energy_arr, L_arr, L0, rtol, atol,
                  _resolve_max_n_slabs(None, integration_method))
              if strategy == 'auto' else None)
''' + s[i + len(old):]
open(p, 'w').write(s); print('ok')
