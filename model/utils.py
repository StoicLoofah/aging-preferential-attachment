from math import *
import gzip
import cPickle as pickle

years = ['2003', '2004', '2005', '2006']
num_bookmarks = {'2003': 106549.0,
                 '2004': 1833989.0,
                 '2005': 16120787.0,
                 '2006': 54627807.0,
                 }
colors = {'2003': 'blue',
          '2004': 'green',
          '2005': 'red',
          '2006': 'cyan'}

#import cPickle as pickle; import gzip; data = pickle.load(gzip.open('', "rb"))

def load(fname):
    file = gzip.open(fname, 'rb')
    data = pickle.load(file)
    file.close()
    return data

def serialize(fname, obj):
    file = gzip.open(fname, 'wb')
    pickle.dump(obj=obj, file=file, protocol=-1)
    file.close()


def tss(p, q):
    max_key = max(max(p.keys()), max(q.keys()))
    diff = 0
    for i in range(max_key):
        v1 = p.get(i, 0)
        v2 = q.get(i, 0)
        diff += pow(v1 - v2, 2)
    return diff

def log_tss(p, q):
    for i in p:
        p[i] = log(p[i])
    for i in q:
        q[i] = log(q[i])
    return tss(p, q)

## {{{ http://code.activestate.com/recipes/66472/ (r1)
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next > end:
            break
        elif inc < 0 and next < end:
            break
        L.append(next)

    return L
## end of http://code.activestate.com/recipes/66472/ }}}
