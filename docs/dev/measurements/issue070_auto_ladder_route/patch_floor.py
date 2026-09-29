p = 'src/magnus/oscprob.py'; s = open(p).read()

def rep(a, b, count=1):
    global s
    c = s.count(a); assert c == count, (c, a[:80]); s = s.replace(a, b)

i0 = s.index('class _PreferLadder:')
i1 = s.index('def _estimated_phase(')
s = s[:i0] + open(__file__.replace('patch_floor.py', 'blk_marker.txt')).read() + s[i1:]
i0 = s.index('def _estimated_phase(')
i1 = s.index('def _osc_prob_hybrid_dispatch(')
s = s[:i0] + open(__file__.replace('patch_floor.py', 'blk_estimate.txt')).read() + s[i1:]

rep('''    if (strategy == 'auto') and _auto_prefers_ladder(H_at_energy, energy_arr, L_arr, L0,
                                                       rtol, atol):
        return _PREFER_LADDER
''', '''    prefer = (_auto_prefers_ladder(H_at_energy, energy_arr, L_arr, L0, rtol, atol)
              if strategy == 'auto' else None)
    if prefer is not None:
        return prefer
''', count=2)
rep('''        prefer_ladder = P_hybrid is _PREFER_LADDER
        if prefer_ladder:
            rtol, atol, scan_kwargs = _with_ladder_margin(rtol, atol, scan_kwargs)
        elif P_hybrid is not NotImplemented:
''', '''        prefer_ladder = isinstance(P_hybrid, _PreferLadder)
        if prefer_ladder:
            rtol, atol, min_n_slabs = P_hybrid.request(rtol, atol, min_n_slabs, max_n_slabs,
                                                       integration_method)
            scan_kwargs = dict(scan_kwargs, rtol=rtol, atol=atol, min_n_slabs=min_n_slabs)
        elif P_hybrid is not NotImplemented:
''', count=2)
rep('''            prefer_ladder = P_scan is _PREFER_LADDER
            if prefer_ladder:
                rtol, atol, scan_kwargs = _with_ladder_margin(rtol, atol, scan_kwargs)
                P_scan = NotImplemented
''', '''            prefer_ladder = isinstance(P_scan, _PreferLadder)
            if prefer_ladder:
                rtol, atol, min_n_slabs = P_scan.request(rtol, atol, min_n_slabs, max_n_slabs,
                                                         integration_method)
                scan_kwargs = dict(scan_kwargs, rtol=rtol, atol=atol, min_n_slabs=min_n_slabs)
                P_scan = NotImplemented
''')
rep('''        # strategy='auto' handed the request to the ladder, which runs at a tenth of the
        # tolerance (issue #70; see osc_prob_matter_std_potential).
        if P_hybrid is _PREFER_LADDER:
            rtol, atol = _with_ladder_margin(rtol, atol)
        elif P_hybrid is not NotImplemented:
''', '''        # strategy='auto' handed the request to the ladder, which runs at a tenth of the
        # tolerance and above a slab floor (issue #70; see osc_prob_matter_std_potential).
        if isinstance(P_hybrid, _PreferLadder):
            rtol, atol, n_floor = P_hybrid.request(rtol, atol, kwargs.get('min_n_slabs'),
                kwargs.get('max_n_slabs'), integration_method)
            kwargs = dict(kwargs, min_n_slabs=n_floor)
        elif P_hybrid is not NotImplemented:
''')
open(p, 'w').write(s)
print('ok')
