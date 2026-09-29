"""Wrap magnus.adiabatic.hybrid_propagator: log every call with atol + rtol < 1e-6."""
import inspect, json, os
import magnus.adiabatic as _ad

def install(context):
    real = _ad.hybrid_propagator
    sig = inspect.signature(real)
    log = os.environ['HP_LOG']
    def wrapped(*args, **kwargs):
        b = sig.bind(*args, **kwargs); b.apply_defaults()
        rtol, atol = b.arguments['rtol'], b.arguments['atol']
        info = b.arguments.get('info')
        own = info is None
        if own:
            b.arguments['info'] = info = {}
        U, windows, cert = real(*b.args, **b.kwargs)
        tight = (rtol or 0.0) + (atol or 0.0) < 1e-6
        if tight:
            with open(log, 'a') as fh:
                fh.write(json.dumps(dict(where=context(), rtol=rtol, atol=atol,
                    n_windows=len(windows), certified=bool(cert),
                    order=b.arguments.get('magnus_exp_order'), pid=os.getpid())) + '\n')
        return U, windows, cert
    _ad.hybrid_propagator = wrapped
    with open(log, 'a') as fh:
        fh.write(json.dumps(dict(hook=_ad.__file__, where=context(), pid=os.getpid())) + '\n')
