import os, sys
sys.path.insert(0, os.environ['HP_DIR'])
import hp_hook
def _ctx():
    ip = get_ipython()
    return '%s | %s' % (os.environ.get('NB_NAME'), ((ip.user_ns.get('_ih') or [''])[-1])[:100].replace('\n', ' / '))
hp_hook.install(_ctx)
