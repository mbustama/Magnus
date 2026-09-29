"""#120 B1 addendum: intermediate margins at the two tight tolerances, same harness and protocol."""
import sys, json
sys.argv = ['x', 'x']
sys.path.insert(0, '/home/mbustamante/Research/magnus/resources/handover_auto_tight')
import tight_measure as tm
tm.TOLS = (1e-9, 1e-12)
tm.CONFIGS = (('H4', 4, None), ('H8', 8, None), ('L4/2', 4, 2.0), ('L4/3', 4, 3.0), ('L8/2', 8, 2.0), ('L8/3', 8, 3.0))
if __name__ == '__main__':
    out = sys.argv_out if hasattr(sys, 'argv_out') else '/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad/tight_b1_margins.jsonl'
    if len(sys.orig_argv) > 2 and sys.orig_argv[-1] == 'score':
        tm.score(out); sys.exit()
    from multiprocessing import Pool
    import verify_b as vb
    jobs = [('scan', n, Es, s) for n, Es, s in tm.scan_cases()] + [('point', n, E, s) for n, E, s in vb.cases()]
    with Pool(4) as pool, open(out, 'w') as fh:
        for rec in pool.imap_unordered(tm.run_case, jobs):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
    print('done', len(jobs))
