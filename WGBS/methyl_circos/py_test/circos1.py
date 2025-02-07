from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset
import numpy as np
np.random.seed(0)

# Load mm10 dataset (https://github.com/moshi4/pycirclize-data/tree/main/eukaryote/mm10)
chr_bed_file = load_eukaryote_example_dataset("mm10")[0]

# Initialize Circos from BED chromosomes
circos = Circos.initialize_from_bed(chr_bed_file, space=2, start=10, end=350, endspace=False)
circos.text("Mus musculus\n(mm10)", size=12)

heatmap_cmaps = ["bwr", "Spectral", "viridis", "plasma", "Reds", "Blues", "Greens", "Greys"]
ColorCycler.set_cmap("gist_rainbow")
colors = ColorCycler.get_color_list(len(circos.sectors))
for sector, color in zip(circos.sectors, colors):
    sector.text(sector.name.replace("chr", ""), size=10)
    # Plot chromosome outer track
    chr_track = sector.add_track((95, 100))
    chr_track.axis(fc=color)
    # Create random test data for heatmap plot
    window_size = 10_000_000
    data_num = int(sector.size // window_size) + 1
    vmin, vmax = 0, 100
    data = np.random.randint(vmin, vmax, data_num)
    # Plot heatmap tracks with various cmap
    track_r_size, r_interval, r_start = 4, 2, 90
    for idx, cmap in enumerate(heatmap_cmaps):
        r_pos = r_start - (track_r_size + r_interval) * idx
        track_r_lim = (r_pos - track_r_size, r_pos)
        track = sector.add_track(track_r_lim)
        track.axis(ec="grey")
        track.heatmap(data, width=window_size, vmin=vmin, vmax=vmax, cmap=cmap)
        # Plot colormap name on center
        if sector.name == circos.sectors[0].name:
            circos.text(cmap, r=track.r_center, size=8)

fig = circos.plotfig()

