from pandas import read_html
from os import mkdir
from os.path import isdir, isfile, join
import requests
from numpy import zeros, empty, object, inf, log
from scipy.interpolate import interp1d
from datetime import datetime
from orthopoly.chebyshev import *

#-------------------------------------------------------------------------------
# INPUT

#url of page with molecule table
molurl = 'https://hitran.org/docs/molec-meta/'
#url of page with isotopologue tables
isourl = 'https://hitran.org/docs/iso-meta/'
#base url for the q/TIPS files
Qurl = 'https://hitran.org/data/Q/'
#directory to store Q files for fast access
Qdir = join('..', 'data', 'Q')
#columns for molecule table
molcols = ['num', 'formula', 'name']
#columns for isotopologue table
isocols = ['gid', 'idx', 'isoform', 'AFGL', 'abundance', 'μ', 'Qref', 'fnQ', 'gi']
#boundaries of temperature range to use for Q interpolation/fitting
Tmin = 25.0
Tmax = 1000.0
#maximum allowable relative error in interpolating function
maxrelerr = 5e-3
#output file
fnout = join('..', 'src', 'molparam.jl')
#log file for info on fitting
fnlog = 'molparam.log'
#TIPS data that must be ignored, specified by the q file name in isourl page
qignore = ['q96.txt']
#show a couple of plots describing the polynomials fits at the end
fitplots = False

#-------------------------------------------------------------------------------
# FUNCTIONS

def isnum(x):

    try:
        float(x)
    except ValueError:
        return(False)
    else:
        return(True)

def parsesci(v):

    if type(v) is str:
        if isnum(v):
            v = float(v)
        else:
            v = v.encode('ascii', 'replace').decode('utf-8')
            n, x = v.split('???')
            v = float(n)*10**float(x[2:])
    return(float(v))

def downloadQ(url):

    r = requests.get(url)
    if r.status_code == 200:
        q = r.content.decode('utf-8').replace('\r', '')
        return(q)
    else:
        return(None)

def parseQ(q, width=4):

    lines = [line.strip() for line in q.splitlines()]
    L = len(lines)
    T, Q = zeros((L,)), zeros((L,))
    for i,line in enumerate(lines):
        T[i] = float(line[:width])
        Q[i] = float(line[width:])
    return(T, Q)

def chebyfit(x, y, Ncheb):

    xmin, xmax = min(x), max(x)
    #first generate a cubic spline to interpolate at chebyshev nodes
    fspline = interp1d(x, y, kind='cubic')
    #setup cheby stuff
    xhat, θ, xc, fc, S = cheby_dct_setup(xmin, xmax, Ncheb)
    #interpolate at cheb nodes
    yc = fspline(xc)
    #compute chev coefficients
    a = fc(yc)
    #evaluate cheb at original points
    ycheb = cheby_sum(x, a, xmin, xmax)
    #compute max abs relative error
    relerr = max(abs(ycheb - y)/abs(y))

    return(relerr, a)

strarray = lambda x: ', '.join([str(i) for i in x])

quotearray = lambda x: ', '.join(['"%s"' % i for i in x])

garray = lambda x: ', '.join(['%.2g' % i for i in x])

write = lambda f, s: f.write(s.encode('utf8'))

#-------------------------------------------------------------------------------
# MAIN

#read the molecule table
dfm = read_html(molurl)[0]
dfm.columns = molcols
#read the isotopologue tables
dfi = read_html(isourl)
for df in dfi:
    df.columns = isocols
#make sure number of molecules matches number of molecule tables
assert len(dfm) == len(dfi), "mismatch between number of molecules"
N = len(dfi)
#parse the strange scientific notation with latin characters
for i in range(N):
    for col in ['abundance', 'Qref']:
        dfi[i][col] = dfi[i][col].apply(parsesci)

#if the q storage directory doesn't exist, make it
if not isdir(Qdir):
    mkdir(Qdir)
#now get all the Q data that's available
for i in range(N):
    df = dfi[i]
    df['T'] = empty((len(df),), dtype=object)
    df['Q'] = empty((len(df),), dtype=object)
    for j in dfi[i].index:
        fnQ = dfi[i].at[j,'fnQ']
        if type(fnQ) is str:
            path = join(Qdir, fnQ)
            if not isfile(path):
                url = join(Qurl, fnQ)
                Q = downloadQ(url)
                print('downloaded: %s' % url)
                print('        to: %s\n' % path)
                with open(path, 'w') as ofile:
                    ofile.write(Q)
            else:
                with open(path, 'r') as ifile:
                    Q = ifile.read()
            T, Q = parseQ(Q)
            df.at[j,'T'] = T
            df.at[j,'Q'] = Q

#fit polnomials to the Q data
print('generating interpolating polynomails for Qref/Q, with T ∈ [%g,%g] K' % (Tmin, Tmax))
f = open(fnlog, 'w')
f.write('molparam.py log\n%s\n\n' % str(datetime.now()))
for i in range(N):
    f.write("molecule %d: %s\n" % (dfm.at[i,'num'], dfm.at[i,'formula']))
    df = dfi[i]
    df['hascheb'] = zeros((len(df),), dtype=bool)
    df['ncheb'] = zeros((len(df),), dtype=int)
    df['coef'] = empty((len(df),), dtype=object)
    df['maxrelerr'] = zeros((len(df),))
    for j in df.index:
        f.write('    isotopologue: %s\n' % df.at[j,'isoform'])
        T, Q = df.at[j,'T'], df.at[j,'Q']
        #ignore this isotopologue if explicitly directed to
        if df.at[j,'fnQ'] in qignore:
            f.write('        ignored at your request!\n')
        #attempt a remedy if there are no T and Q data
        elif T is None and Q is None:
            #if no q data, try to copy from most abundant isotopologue
            if df.at[0,'ncheb'] != 0:
                df.at[j,'ncheb'] = df.at[0,'ncheb']
                df.at[j,'coef'] = df.at[0,'coef']
                df.at[j,'maxrelerr'] = df.at[0,'maxrelerr']
                df.at[j,'hascheb'] = True
                f.write('        no Q data, copied fit from %s\n' % df.at[0,'isoform'])
            #otherwise, it's a bust
            else:
                f.write('        no Q data, cannot generate a fit for Qref/Q\n')
        #data are there, construct the fit!
        else:
            T, Q = df.at[j,'T'], df.at[j,'Q']
            Qref = df.at[j,'Qref']
            #restrict the temperature range and define the curve to fit
            m = (T >= Tmin) & (T <= Tmax)
            x = T[m]
            y = Q[m]/Qref
            #make increasingly accurate cheby fits until its good enough
            n = 3 #require at least 3 terms
            e = inf
            while e > maxrelerr:
                e, a = chebyfit(x, y, n)
                n += 1
            n -= 1
            #store in the frame
            df.at[j,'ncheb'] = n
            df.at[j,'coef'] = a
            df.at[j,'hascheb'] = True
            df.at[j,'maxrelerr'] = e
            f.write('        max relative err: %g\n' % e)
            f.write('        polynomial terms: %d\n' % n)

    f.write('\n')
f.close()
print('file written: %s' % fnlog)

#print a note about the max number of coefficients and largest error
print('The largest number of coefficients is %d' % max(df['ncheb'].max() for df in dfi))
print('The largest relative error anywhere is %g' % max(df['maxrelerr'].max() for df in dfi))

#write to Julia!
with open(fnout, 'wb') as f:

    #write interpolation boundaries as constants
    write(f, 'const TMIN = %f\n' % Tmin)
    write(f, 'const TMAX = %f\n' % Tmax)

    #write molparam structs
    write(f, 'const MOLPARAM = MolParam[\n')
    for i in range(1, dfm['num'].max() + 1):
        #the molecule number must be present
        if i in dfm['num'].values:
            #get the index of the molecule number
            idx = dfm['num'][dfm['num'] == i].index[0]
            #get the proper isotopologue table
            df = dfi[idx]
            #write all the info
            write(f, '  MolParam(\n')
            write(f, '    #1, molecule number\n')
            write(f, '    %d,\n' % i)
            write(f, '    #2, molecule formula\n')
            write(f, '    "%s",\n' % dfm.at[idx,'formula'])
            write(f, '    #3, molecule name\n')
            write(f, '    "%s",\n' % dfm.at[idx,'name'])
            write(f, '    #4, global isotopologue numbers\n')
            write(f, '    Int64[%s],\n' % strarray(df['gid']))
            write(f, '    #5, isotopologue formulae\n')
            write(f, '    String[%s],\n' % quotearray(df['isoform']))
            write(f, '    #6, AFGL code\n')
            write(f, '    Int64[%s],\n' % strarray(df['AFGL']))
            write(f, '    #7, abundance fractions\n')
            write(f, '    Float64[%s],\n' % strarray(df['abundance']))
            write(f, '    #8, molecular masses [kg/mole]\n')
            write(f, '    Float64[%s],\n' % strarray(df['μ']/1e3))
            write(f, '    #9, Qref\n')
            write(f, '    Float64[%s],\n' % strarray(df['Qref']))
            write(f, '    #10, has interpolating chebyshev expansion?\n')
            write(f, '    Bool' + ('[%s],\n' % strarray(df['hascheb'])).lower())
            write(f, '    #11, lengths of interpolating chebyshev expansion\n')
            write(f, '    Int64[%s],\n' % strarray(df['ncheb']))
            write(f, '    #12, maximum relative errors of interpolation\n')
            write(f, '    Float64[%s],\n' % garray(df['maxrelerr']))
            write(f, '    #13, chebyshev expansion coefficients\n')
            write(f, '    Vector{Float64}[\n')
            for j,coef in enumerate(df['coef']):
                if coef is not None:
                    write(f, '      Float64[%s]' % strarray(coef))
                else:
                    write(f, '      Float64[]')
                if j == len(df['coef']) - 1:
                    write(f, '\n    ]\n')
                else:
                    write(f, ',\n')
            #close the struct definition
            if i == dfm['num'].max():
                write(f, '  )\n')
            else:
                write(f, '  ),\n')
        else:
            #no molecule with that number
            write(f, '  MolParam(), #no molecule has been assigned the number %d\n' % i)
    #close the outer array
    write(f, ']')
print('file written: %s\n' % fnout)

if fitplots:
    from pandas import concat
    df = concat(dfi, axis=0)

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2,1)
    axs[0].hist(df['maxrelerr'], color='gray')
    axs[0].set_xlabel('Maximum Relative Error')
    axs[0].set_ylabel('Counts')
    axs[1].hist(df['ncheb'][df['ncheb'] > 0],
        color='gray',
        range=(0.5, df['ncheb'].max() + 0.5),
        bins=df['ncheb'].max())
    axs[1].set_xlabel('Chebyshev Terms')
    axs[1].set_ylabel('Counts')
    fig.tight_layout()
    plt.show()
