# -*- coding: utf-8 -*-

def filter_chan(channel):
    # you can check the channel name here and only return True if it's in some whitelist, etc. to build
    # a model which consists of some sub-channels only. Here, we keep everything:
    return True


# a "signal process group" is a list of process names to be scaled *together* by the signal strength parameter "beta_signal"
signal_process_groups = {'ttZttWZ': ['ttZ', 'ttWZ']}
model = higgs_datacard.build_model('datacard_7channels.txt', filter_channel = filter_chan)
model.signal_process_groups = signal_process_groups


# calculate profile likelihood interval on data:
res = pl_interval(model, 'data', 1)['ttZttWZ']
#print res
print "profile likelihood interval of the signal strength parameter: %.3f   +%.3f  -%.3f" % (res[0.0][0], -res[0.0][0] + res[cl_1sigma][0][1], res[0.0][0] - res[cl_1sigma][0][0])

# calculate the approximate significance (Z-value):
zvalues = zvalue_approx(model, 'data', n = 1)['ttZttWZ']
print "approximate Z-value: %.3f" % zvalues[0]


res = pl_interval(model, 'data', 1)['ttZttWZ']
print "complete pl interval for beta_signal: %.3f   +%.3f  -%.3f" % (res[0.0][0], -res[0.0][0] + res[cl_1sigma][0][1], res[0.0][0] - res[cl_1sigma][0][0])


res = mle(model, 'data', 1)['ttZttWZ']
values = {}
for p in model.get_parameters(signal_process_groups['ttZttWZ']): values[p] = res[p][0][0]
fix_nuisance_pars = get_fixed_dist_at_values(values)
res = pl_interval(model, 'data', 1, nuisance_constraint = fix_nuisance_pars)['ttZttWZ']
print "pl interval for beta_signal, fixing nuisances: %.3f   +%.3f -%.3f" % (res[0.0][0], -res[0.0][0] + res[cl_1sigma][0][1], res[0.0][0] - res[cl_1sigma][0][0])


# calculate the significance using toys (will take a long time!):
#discovery(model)


model_summary(model)
report.write_html('htmlout')

