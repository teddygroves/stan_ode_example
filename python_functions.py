import pandas as pd
import arviz
import numpy as np
import pystan
import pickle
from hashlib import md5
from matplotlib import pyplot as plt


def ppc_plot(infd: arviz.InferenceData, name=None, sharey=False, labels=False, ylim=None):
    obs = infd.observed_data.to_dataframe()
    posterior_pred = infd.posterior_predictive['metabolite_pred'].to_dataframe()
    prior_pred = infd.prior_predictive['metabolite_pred'].to_dataframe()

    g_posterior = posterior_pred.groupby(['experiments', 'metabolites'])
    g_prior =  prior_pred.groupby(['experiments', 'metabolites'])

    f, axes = plt.subplots(1, 2, figsize=[15, 5], sharey=sharey)

    for ax, g, name in zip(axes, [g_prior, g_posterior], ['prior', 'posterior']):
        lower = g['metabolite_pred'].quantile(0.1).rename('lower')
        mean = g['metabolite_pred'].mean().rename('mean')
        upper = g['metabolite_pred'].quantile(0.9).rename('upper')

        xx = np.linspace(mean.min(), mean.max(), 10)
        scatter = ax.scatter(mean.reindex(obs.index), obs, marker='x')
        vlines = ax.vlines(mean, lower, upper, color='tab:orange', zorder=0, label='10%-90% credible interval')
        y_equals_x_line = ax.plot(xx, xx, color='r', linestyle='--', label='y=x')
        if ylim is not None:
            ax.set_ylim(ylim)
        if labels:
            for i, m in g[['metabolite_pred']].mean().iterrows():
                ax.text(m, obs.loc[i], '-'.join(i), fontsize=8, horizontalalignment='center')
        text = ax.set(xlabel=f'{name.capitalize()} mean concentration',
                      ylabel='Observed concentration')
        title = ax.set_title(f'{name.capitalize()} predictive check', y=0.8)
    plt.figlegend([vlines, y_equals_x_line[0], scatter], 
                  ['10%-90% credible interval', 'y=x', 'Metabolite/experiment'],
                  ncol=3, loc='lower center', bbox_to_anchor=[0.45, 0.9])
    if name is None:
        name = 'checks'
    plt.savefig(f'output/{name}.png')


def StanModel_cache(file, model_name=None, **kwargs):
    with open(file, 'r') as f:
        model_code = f.read()
    code_hash = md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached_stan_models/{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached_stan_models/{}-{}.pkl'.format(model_name, code_hash)
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code)
        with open(cache_fn, 'wb') as f:
            pickle.dump(sm, f)
    else:
        print("Using cached StanModel")
    return sm
