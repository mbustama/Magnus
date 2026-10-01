import warnings, numpy as np, magnus.oscprob as oscprob, magnus.globaldefs as gd
for c,d in [(-0.272,20.0),(-0.552,1000.0),(-0.872,2800.0)]:
  with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    P = oscprob.osc_prob_3nu_earth(2.5*gd.UNIT_MEV, costhz=c, source_depth=d*gd.UNIT_KM,
        detector_depth=1.4*gd.UNIT_KM, nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, density_matter_ocean=2.65, average=True)
  print(c, P, [str(x.message)[str(x.message).find('largest'):][:70] for x in w if 'largest' in str(x.message)])
