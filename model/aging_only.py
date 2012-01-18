import networkx as nx
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import scipy.stats as stats
import random
from math import *
import cPickle as pickle
from datetime import *
import gzip
import itertools
from utils import *

gamma = 0
alpha = 0


def parse_starts_and_dates(years=years):
    for year in years:
        data = load(year + 'cascades.pkl.gz')
        starts = {}
        dates = {}
        first_date = datetime(int(year), 1, 1).toordinal()
        for all_dates in data.itervalues():
            cur_start_date = all_dates[0].toordinal() - first_date
            if cur_start_date in starts:
                starts[cur_start_date] += 1
            else:
                starts[cur_start_date] = 1
            for date in all_dates[1:]:
                cur_date = date.toordinal() - first_date
                if cur_date in dates:
                    dates[cur_date] += 1
                else:
                    dates[cur_date] = 1
        starts = sorted(starts.items(), key=lambda d: d[0])
        dates = sorted(dates.items(), key=lambda d: d[0])
        data = {'starts': starts, 'dates': dates}
        serialize(year + 'starts_and_dates.pkl.gz', data)
        print('Done with ' + year)

def generate_ao_network(year, alpha=alpha,
                           to_save=False, compare=False):
    alpha = -alpha
    data = load('starts_and_dates/{0}starts_and_dates.pkl.gz'.format(year))

    start_dates = []
    for date, size in data['starts']:
        start_dates.extend([date for i in range(size)])
    num_start_dates = len(start_dates)
    sizes = [1] * num_start_dates
    cur_limit = 0

    delay_dist = {}

    for i, pair in enumerate(data['dates']):
        after_date = pair[0]
        while cur_limit < num_start_dates and \
                start_dates[cur_limit] <= after_date:
            cur_limit += 1
        factors = [0] * cur_limit
        factor_sum = 0
        for j in range(cur_limit):
            factor_sum += pow(after_date - start_dates[j] + 1, alpha)
            factors[j] = factor_sum
        for j in range(pair[1]):
            random_value = random.random() * factor_sum
            if random_value < factors[0]:
                sizes[0] += 1
            else:
                left = 0
                right = cur_limit
                while left < right:
                    mid = (left + right) / 2
                    if factors[mid] < random_value:
                        left = mid + 1
                    elif factors[mid - 1] > random_value:
                        right = mid - 1
                    else:
                        break
                random_index = (left + right) / 2
                sizes[random_index] += 1
                delay = after_date - start_dates[random_index]
                if delay in delay_dist:
                    delay_dist[delay] += 1
                else:
                    delay_dist[delay] = 1
    deg_dist = {}
    for s in sizes:
        if s in deg_dist:
            deg_dist[s] += 1
        else:
            deg_dist[s] = 1

    vals = sorted(deg_dist.items(), key=lambda item: item[0])
    deg = [val[0] for val in vals]
    count = [val[1] for val in vals]

    estimate = False
    if estimate:
        alpha_mle = 0
        num_values = 0
        x_min = 1
        for key, value in zip(deg, count):
            alpha_mle += value * log(key / x_min)
            num_values += value
        alpha_mle = 1 + num_values * pow(alpha_mle, -1)
        error = (alpha_mle - 1)/(sqrt(num_values))
        print('mle estimate of alpha: {0}\n'.format(alpha_mle))
        print('mle standard error: {0}\n\n'.format(error))

    if to_save:
        serialize('deg_dist/{0}_ao_{1}.pkl.gz'.format(year, alpha), deg_dist)
        serialize('delay_dist/{0}_ao_{1}.pkl.gz'.format(year, alpha), delay_dist)

    if compare:
        actual_deg_dist = load("deg_dist/{0}.pkl.gz".format(year))
        actual_delay_dist = load("delay_dist/{0}.pkl.gz".format(year))
        deg_tss = log_tss(actual_deg_dist, deg_dist)
        delay_tss = log_tss(actual_delay_dist, delay_dist)
        return deg_tss, delay_tss

    '''
def plot_aging_nets():
    year = '2003'
    gamma = 1.45
    alpha = -.5

    deg_dist = load('deg_dist/{0}_model_g{1}_a{2}.pkl.gz'.format(
            year, gamma, alpha))
    actual_deg_dist = load("deg_dist/{0}.pkl.gz".format(year))

    plt.clf()

    vals = sorted(deg_dist.items(), key=lambda item: item[0])
    deg = [val[0] for val in vals]
    count = [(val[1] / num_bookmarks[year]) for val in vals]
    plt.loglog(deg, count, label='model')

    vals = sorted(actual_deg_dist.items(), key=lambda item: item[0])
    deg = [val[0] for val in vals]
    count = [(val[1] / num_bookmarks[year]) for val in vals]
    plt.loglog(deg, count, label='actual')

    plt.xlabel('degree')
    plt.ylabel('P(k)')
    plt.legend()
    plt.savefig('plots/aging_combined_{0}.png'.format(year))

    delay_dist = load('delay_dist/{0}_model_g{1}_a{2}.pkl.gz'.format(
            year, gamma, alpha))
    plt.clf()
    items = sorted(delay_dist.items(), key=lambda d: d[0])
    delay = [item[0] for item in items]
    count = [(item[1] / num_bookmarks[year]) for item in items]
    plt.loglog(delay, count, label='model')

    actual_delay_dist = load("delay_dist/{0}.pkl.gz".format(year))
    items = sorted(actual_delay_dist.items(), key=lambda d: d[0])
    delay = [item[0] for item in items]
    count = [(item[1] / num_bookmarks[year]) for item in items]
    plt.loglog(delay, count, label='actual')

    plt.xlabel('age (days)')
    plt.ylabel('P(k)')
    plt.legend()
    plt.savefig('plots/aging_combined_delay_{0}.png'.format(year))

'''
def search_for_params(year, min_alpha, max_alpha, inc_alpha):
    fout = open('results/{0}ao_results{1}.csv'.format(year, min_alpha), 'w')

    alphas = frange(min_alpha, max_alpha, inc_alpha)

    for alpha in alphas:
        print(datetime.now())
        deg_tss, delay_tss = generate_ao_network(
            year, alpha, to_save=True, compare=True)
        print('alpha={1}\ndegree TSS={2}\ndelay TSS={3}\ncombined TSS={4}\n'.format(0, alpha, deg_tss, delay_tss, deg_tss + delay_tss))
        fout.write('{1},{2},{3},{4}\n'.format(
                0, alpha, deg_tss, delay_tss, deg_tss + delay_tss))

    fout.close()


def test_inputs():
    while True:
        cur_alpha = input('Alpha: ')
        deg_tss, delay_tss = generate_ao_network(
            '2003', cur_alpha, plot=True, compare=True)

#parse_starts_and_dates()
#test_inputs()
#search_for_params()
#plot_aging_nets()
search_for_params('2004', 0, 1, .05)

