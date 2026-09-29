import time; t0 = time.perf_counter()
import sys, json
import numpy as np, magnus.globaldefs as gd, magnus.oscprob as oscprob
t1 = time.perf_counter()
exec(compile(open(sys.argv[1]).read(), 'listing', 'exec'), {})
t2 = time.perf_counter()
print(json.dumps(dict(imports=t1 - t0, calls=t2 - t1)))
