"""pytest plugin: the same hook, with the running test's node id as context."""
import os
CURRENT = {'id': 'collect'}
def pytest_configure(config):
    import hp_hook
    hp_hook.install(lambda: CURRENT['id'])
def pytest_runtest_setup(item):
    CURRENT['id'] = item.nodeid
