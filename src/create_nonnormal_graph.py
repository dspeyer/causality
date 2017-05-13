import datareader
from scipy.stats.mstats import zscore
import utils

d = datareader.Data()
d.read_from_files()
pl = d.get_pid_list(lambda(p): ('nod2' in p and 
                                   p.health in ['control', 'ileal CD'] and
                                   len(p.fractions)>0))
data=d
ent={}
for species in data.bacteria:
    vals = data.get_data(pl, 'fractions', species, bucketizer=lambda(x):int(log(x,10)))
    ent[species] = utils.entropy(vals)

spe=keys(ent)
spe=ent.keys()
spe.sort(key=lambda x: ent[x], reverse=True)

