import os, json
try:
    import magnus.oscprob as _op
    _orig = _op._auto_prefers_ladder
    def _logged(*a, **k):
        r = _orig(*a, **k)
        if r is not None:
            ip = get_ipython()
            with open(os.environ['ROUTE_LOG'], 'a') as fh:
                fh.write(json.dumps(dict(nb=os.environ.get('NB_NAME'), exec_count=ip.execution_count,
                                         cell=(ip.user_ns.get('_ih') or [''])[-1][:120],
                                         floor=r.min_n_slabs)) + '\n')
        return r
    _op._auto_prefers_ladder = _logged
    with open(os.environ['ROUTE_LOG'], 'a') as fh:
        fh.write(json.dumps(dict(nb=os.environ.get('NB_NAME'), hook=_op.__file__)) + '\n')
except Exception as exc:
    print('route hook FAILED', exc)
