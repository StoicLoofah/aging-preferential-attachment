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

def generate_aging_network(year, gamma=gamma, alpha=alpha,
                           plot=False, compare=False):
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
            factors[j] = pow(sizes[j], gamma) * \
                pow(after_date - start_dates[j] + 1, alpha)
            factor_sum += factors[j]
        for j in range(pair[1]):
            random_value = random.random() * factor_sum
            random_index = cur_limit - 1
            cur_sum = 0
            for i, e in enumerate(factors):
                cur_sum += e
                if random_value < cur_sum:
                    random_index = i
                    break
            sizes[random_index] += 1
            delay = after_date - start_dates[random_index]
            if delay in delay_dist:
                delay_dist[delay] += 1
            else:
                delay_dist[delay] = 1
            new_factor = pow(sizes[random_index], gamma) * \
                pow(delay + 1, alpha)
            factor_sum += new_factor - factors[random_index]
            factors[random_index] = new_factor
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

    if plot:
        serialize('deg_dist/{0}_model_g{1}_a{2}.pkl.gz'.format(
                year, gamma, alpha), deg_dist)
        plt.clf()
        plt.loglog(deg, count, label='model')

        actual_deg_dist = load("deg_dist/{0}.pkl.gz".format(year))
        vals = sorted(actual_deg_dist.items(), key=lambda item: item[0])
        deg = [val[0] for val in vals]
        count = [val[1] for val in vals]
        plt.loglog(deg, count, label='actual')

        plt.xlabel('size')
        plt.ylabel('count')
        plt.legend()
        plt.savefig('aging_combined_{0}_g{1}_a{2}_.png'.format(
                year, gamma, alpha))

        serialize('delay_dist/{0}_model_g{1}_a{2}.pkl.gz'.format(
                year, gamma, alpha), delay_dist)
        plt.clf()
        items = sorted(delay_dist.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        degree = [item[1] for item in items]
        plt.loglog(delay, degree, label='model')

        actual_delay_dist = load("delay_dist/{0}.pkl.gz".format(year))
        items = sorted(actual_delay_dist.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        degree = [item[1] for item in items]
        plt.loglog(delay, degree, label='actual')

        plt.xlabel('age')
        plt.ylabel('count')
        plt.legend()
        plt.savefig('aging_combined_delay{0}_g{1}_a{2}.png'.format(
                year, gamma, alpha))

    if compare:
        deg_tss = log_tss(actual_deg_dist, deg_dist)
        delay_tss = log_tss(actual_delay_dist, delay_dist)
        return deg_tss, delay_tss

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


def search_for_params(year, min_gamma, max_gamma, inc_gamma, min_alpha, max_alpha, inc_alpha):
    fout = open('results/{0}aging_results{1}{2}.csv'.format(year, min_gamma, min_alpha), 'w')

    gammas = frange(min_gamma, max_gamma, inc_gamma)
    alphas = frange(min_alpha, max_alpha, inc_alpha)

    for alpha in alphas:
        for gamma in gammas:
            print(datetime.now())
            deg_tss, delay_tss = generate_aging_network(
                year, gamma, alpha, plot=True, compare=True)
            print('gamma={0}\nalpha={1}\ndegree TSS={2}\ndelay TSS={3}\ncombined TSS={4}\n'.format(gamma, alpha, deg_tss, delay_tss, deg_tss + delay_tss))
            fout.write('{0},{1},{2},{3},{4}\n'.format(
                    gamma, alpha, deg_tss, delay_tss, deg_tss + delay_tss))

    fout.close()


def test_inputs():
    while True:
        cur_gamma = input('Gamma: ')
        cur_alpha = input('Alpha: ')
        deg_tss, delay_tss = generate_aging_network(
            '2003', cur_gamma, cur_alpha, plot=True, compare=True)

#parse_starts_and_dates()
#test_inputs()
#search_for_params()
plot_aging_nets()
