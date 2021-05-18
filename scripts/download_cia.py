from pandas import read_html
from os import mkdir
from os.path import join, isdir
import requests

#-------------------------------------------------------------------------------
# INPUT

#url of Hitran page with CIA table and download links
ciaurl = 'https://hitran.org/cia/'
#base url for downloading data
dataurl = 'https://hitran.org/data/CIA'
#output directory
dirout = join('..', 'data', 'cia')

#-------------------------------------------------------------------------------
# FUNCTIONS

#create the output directory if need be
if not isdir(dirout):
    mkdir(dirout)
print('download directory: %s' % dirout)

#get the table of available data
df = read_html(ciaurl)[0]

#download each file
for fn in df['Filename']:
    url = '/'.join([dataurl, fn])
    r = requests.get(url)
    if r.status_code == 200:
        path = join(dirout, fn)
        with open(path, 'wb') as f:
            f.write(r.content)
            print('+ %s' % fn)
    else:
        print('request failed for url: %s' % url)
