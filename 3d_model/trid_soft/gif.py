#!/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3
import sys
import glob
from PIL import Image

infile = sys.argv[1]
outfile = sys.argv[2]
# Create the frames
#imgs = sorted(glob.glob("fig_*.png"))
imgs = sorted([i.strip() for i in open(infile).readlines()])
frames = []
for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 
# Save into a GIF file that loops forever
frames[0].save(outfile, format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=60, loop=0)
