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

def serialize_cascades():
    for year in years:
        fin = gzip.open(year + 'cascades.csv.gz', 'r')
        resources = {}
        for line in fin:
            data = line[0:-1].split(',')
            resources[int(data[0])] = sorted(map(
                lambda date_string: datetime.strptime(date_string, '%Y-%m-%d'),
                data[1:]))
        fin.close()

        serialize('cascades/' + year + '.pkl.gz', resources)
        print('Done with ' + year)


def find_cascades():
    '''
    turns raw data into cascades
    '''
    for year in years:
        fin = open(year + '.dat', 'r')
        resources = {}
        for line in fin:
            data = line.split('\t')
            if data[2] in resources:
                resources[data[2]].append(data[0])
            else:
                resources[data[2]] = [data[0]]
        fin.close()

        fout = open(year + 'cascades.csv', 'w')
        for resource, vals in resources.iteritems():
            output = (resource + ',' + (','.join(vals))).replace('\n', '')
            fout.write(output + '\n')
        fout.close()

        d = degree_distribution(resources)
        degree_dists = sorted(d.items(), lambda x,y: x[0] < y[0])

        fout = open(year + 'degdist.csv', 'w')
        for key, value in degree_dists:
            fout.write(str(key) + ',' + str(value) + '\n')
        fout.close()
        print('Done with ' + year)


def degree_distribution(resources):
    d = {}
    for r in resources:
        cur_d = len(resources[r])
        if cur_d in d:
            d[cur_d] += 1
        else:
            d[cur_d] = 1
    return d


def extract_degree_dists(years=years):
    '''
    turned old csvs into seralized objects
    '''
    for year in years:
        deg_dist = {}
        fin = open('results/' + year + 'degdist.csv', 'r')
        for line in fin:
            comma = line.find(',')
            deg_dist[int(line[0:comma])] = int(line[comma+1:])
        serialize('deg_dist/{0}.pkl.gz'.format(year), deg_dist)


def plot_degree_dists(years=years):
    '''
    plots degree distribution and estimates params
    '''
    estimate = False
    plt.clf()
    for year in years:
        data = load('deg_dist/{0}.pkl.gz'.format(year))
        items = sorted(data.items(), key=lambda d: d[0])
        deg = [val[0] for val in items]
        count = [(val[1] / num_bookmarks[year]) for val in items]
        plt.loglog(deg, count, label=year)

        if estimate:
            log_deg = map(log, deg)
            log_count = map(log, count)
            slope, intercept, r_value, p_value, std_err = \
                stats.linregress(log_deg, log_count)
            fout.write('slope: {0}\nintercept: {1}\nr-value: {2}\np-value: {3}\nstd err: {4}\n'.format(slope, intercept, r_value, p_value, std_err))

            alpha_mle = 0
            num_values = 0
            x_min = 1
            for key, value in zip(deg, count):
                alpha_mle += value * log(key / x_min)
                num_values += value
            alpha_mle = 1 + num_values * pow(alpha_mle, -1)
            error = (alpha_mle - 1)/(sqrt(num_values))
            #fout.write('mle estimate of alpha: {0}\n'.format(alpha_mle))
            #fout.write('mle standard error: {0}\n\n'.format(error))
    #fout.close()
    plt.legend()
    plt.xlabel('degree')
    plt.ylabel('P(k)')
    plt.savefig('plots/degree_dist_p.png')


def get_delay(years=years, suffix=''):
    '''
    this function looks at the cascades and figures out the distribution
    over the time it takes for the followups to come. Not entirely sure
    what this will turn out to be
    '''

    fout = open('results/delay_average{0}.csv'.format(suffix), 'w')
    for year in years:
        data = load('cascades/{0}{1}.pkl.gz'.format(year, suffix))

        count = 0
        total = 0
        delay_dist = {}
        for resource, vals in data.iteritems():
            size = len(vals)
            start_date = vals[0]
            for value in vals[1:]:
                diff = (value - start_date).days
                if diff in delay_dist:
                    delay_dist[diff] += 1
                else:
                    delay_dist[diff] = 1

                total += diff
                count += 1

        fout.write('{0},{1}\n'.format(year, total / count))

        serialize('delay_dist/{0}{1}.pkl.gz'.format(year, suffix), delay_dist)
        print('Done with ' + year)
    fout.close()

def plot_delay(years=years, suffix=""):
    '''
    takes the result from the get_delay function and turns it into a plot
    '''

    plt.clf()
    for year in years:
        data = load('delay_dist/{0}{1}.pkl.gz'.format(year, suffix))
        items = sorted(data.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        degree = [(item[1] / num_bookmarks[year]) for item in items]
        plt.loglog(delay, degree, label=year)
    plt.legend(loc=3)
    plt.xlabel('age (days)')
    plt.ylabel('P(k)')
    plt.savefig('plots/delay_plot_p{0}.png'.format(suffix))

def plot_delay_combined(years=years):
    '''
    takes the result from the get_delay function and turns it into a plot
    '''

    plt.clf()
    suffix=''
    for year in years:
        data = load('delay_dist/{0}{1}.pkl.gz'.format(year, suffix))
        items = sorted(data.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        degree = [(item[1] / num_bookmarks[year]) for item in items]
        plt.loglog(delay, degree, label=year)
    suffix='_random'
    for year in years:
        data = load('delay_dist/{0}{1}.pkl.gz'.format(year, suffix))
        items = sorted(data.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        degree = [(item[1] / num_bookmarks[year]) for item in items]
        plt.loglog(delay, degree, label=(year + suffix))
    plt.legend()
    plt.xlabel('age (days)')
    plt.ylabel('P(k)')
    plt.savefig('plots/delay_plot_p_combined.png')


def find_time(suffix=''):
    fout = open('results/date_average{0}.csv'.format(suffix), 'w')
    for year in years:
        fname = 'cascades/{0}{1}.pkl.gz'.format(year, suffix)
        data = pickle.load(gzip.open(fname, "rb"))

        count = 0
        total = 0
        date_dist = {}
        for resource, vals in data.iteritems():
            size = len(vals)
            if size > 1:
                average = sum([(val - vals[0]).days for val in vals[1:]]) / \
                    float(size - 1)
                if size in date_dist:
                    date_dist[size].append(average)
                else:
                    date_dist[size] = [average]
                total += average
                count += 1
        fout.write('{0},{1}\n'.format(year, total / count))

        serialize('date_dist/{0}{1}.pkl.gz'.format(year, suffix), date_dist)
        fname = year + 'date_dist{0}.pkl.gz'.format(suffix)
        pickle.dump(obj=date_dist, file=gzip.open(fname, 'wb'), protocol=-1)
        print('Done with ' + year)
    fout.close()


def get_time_plots(suffix=''):
    fout = open('results/time_distribution{0}.txt'.format(suffix), 'w')
    plt.clf()
    for year in years:
        date_dist = load('date_dist/{0}{1}.pkl.gz'.format(year, suffix))
        vals = sorted(date_dist.items(), key=lambda d: d[0])
        deg = [log(val[0]) for val in vals]
        average = [(sum(val[1]) / float(len(val[1]))) for val in vals]
        plt.scatter(deg, average, label=year, c=colors[year], marker='o')

        total = 0.0
        count = 0
        for key, value in vals:
            count += len(value)
            total += sum(value)

        fout.write(year + '\n')
        fout.write('Average time: {0}\n'.format(total / count))
        print('Done with ' + year)
    fout.close()
    plt.legend()
    plt.xlabel('degree (log)')
    plt.ylabel('average age (days)')
    plt.savefig('plots/time_plot{0}.png'.format(suffix))


def filter_cascades():
    for year in years:
        fname = year + 'cascades.pkl.gz'
        data = pickle.load(gzip.open(fname, "rb"))
        for resource in data.keys():
            if len(data[resource]) == 1:
                del data[resource]
        fname = year + 'cascades_filtered.pkl.gz'
        pickle.dump(obj=data, file=gzip.open(fname, 'wb'), protocol=-1)


def create_random_networks():
    '''
    shuffles everything, then puts them back in as cascades
    '''
    for year in years:
        data = load('cascades_filtered/' + year + 'cascades_filtered.pkl.gz')
        dates = list(itertools.chain.from_iterable(data.itervalues()))
        random.shuffle(dates)
        i = 0
        for resource in data.iterkeys():
            new_i = i + len(data[resource])
            data[resource] = sorted(dates[i:new_i])
            i = new_i
        serialize('cascades/' + year + '_random.pkl.gz', data)
        print('Done with ' + year)

'''
def sort_random_cascades():
    for year in years:
        fname = year + 'cascades_random.pkl.gz'
        data = pickle.load(gzip.open(fname, "rb"))
        for val in data.itervalues():
            val.sort()
        pickle.dump(obj=data, file=gzip.open(fname, 'wb'), protocol=-1)
        print('Done with ' + year)
'''

def plot_models():
    '''
    year = '2003'
    lines = [
        ['2003.pkl.gz', 'actual',],
        ['2003_pa_1.0.pkl.gz', 'pref only', ],
        ['2003_ao_-0.25.pkl.gz', 'aging only', ],
        ['2003_model_g1.45_a-0.5.pkl.gz', 'combined', ],
        ]
        '''
    year = '2004'
    lines = [
        ['2004.pkl.gz', 'actual',],
        ['2004_pa_1.0.pkl.gz', 'pref only', ],
        ['2004_ao_0.pkl.gz', 'aging only', ],
        #['2004_model_g1.45_a-0.5.pkl.gz', 'combined', ],
        ]
    plt.clf()
    for line in lines:
        deg_dist = load('deg_dist/' + line[0])
        vals = sorted(deg_dist.items(), key=lambda item: item[0])
        deg = [val[0] for val in vals]
        count = [(val[1] / num_bookmarks[year]) for val in vals]
        plt.loglog(deg, count, label=line[1])

    plt.xlabel('degree')
    plt.ylabel('P(k)')
    plt.legend()
    plt.savefig('plots/aging_combined_{0}.png'.format(year))

    plt.clf()
    for line in lines:
        delay_dist = load('delay_dist/' + line[0])
        items = sorted(delay_dist.items(), key=lambda d: d[0])
        delay = [item[0] for item in items]
        count = [(item[1] / num_bookmarks[year]) for item in items]
        plt.loglog(delay, count, label=line[1])

    plt.xlabel('age (days)')
    plt.ylabel('P(k)')
    plt.legend()
    plt.savefig('plots/aging_combined_delay_{0}.png'.format(year))



#plot_degree_dists()
#find_cascades()
#serialize_cascades()
#find_time()
#filter_cascades()
#create_random_networks()
#sort_random_cascades()
#find_time('_random')
#get_time_plots()
#get_time_plots('_random')
#get_delay()
#get_delay(years, '_random')
#plot_delay(years)
#plot_delay(years, '_random')
#plot_delay_combined(years)
#extract_degree_dists()
plot_models()
