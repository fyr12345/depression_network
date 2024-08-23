
import brainsmash

#whole-brain
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.eval import sampled_fit
coord_file = r"D:\shiyan\lingmo\sample_in_mask_coord.txt"
output_dir = r"D:\shiyan\lingmo"
brain_map = r"D:\shiyan\lingmo\t4.txt"
filenames = volume(coord_file, output_dir)  #to generate memory-mapped  diatance matrix
kwargs = {'ns':100, 'knn': 900, 'pv': 40, 'nh': 12, 'kernel': 'uniform'}
sampled_fit(brain_map, filenames['D'], filenames['index'], nsurr=10, **kwargs) #visually inspect the variogram fit
from brainsmash.mapgen.sampled import Sampled
gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], **kwargs) #to randomly generate surrogate brain maps with SA that is matched to our target brain map
surrogate_maps = gen(n=10000)
import numpy as np
np.savetxt(r"D:\shiyan\lingmo\t4.csv", surrogate_maps, delimiter="," )

