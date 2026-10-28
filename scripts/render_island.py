"""3D frames and a GIF of an island run with PyVista (pip install pyvista),
for machines without ParaView. Same scene as paraview_island.py, cropped to
the 18 m around the island. Heights exaggerated 10x.

    python scripts/render_island.py results/island_c/island.pvd figures/island_c.gif
"""
import os
import sys
import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
import vtu

pvd = sys.argv[1] if len(sys.argv) > 1 else "results/island_c/island.pvd"
gif = sys.argv[2] if len(sys.argv) > 2 else "figures/island_c.gif"
S = 10.0
pv.OFF_SCREEN = True
FOAM = LinearSegmentedColormap.from_list("foam", [(0.16, 0.47, 0.84), (0.92, 0.96, 1.0)])

frames = []
for t, path in vtu.series(pvd):
    m = pv.read(path).cell_data_to_point_data(pass_cell_data=True).clip_box([4, 22, 4.8, 22.8, -10, 10], invert=False)
    m["speed"] = np.hypot(m["u"], m["v"])
    ground = m.warp_by_scalar("b", factor=S)
    water = m.threshold(1e-3, scalars="h", preference="cell").warp_by_scalar("eta", factor=S)

    p = pv.Plotter(off_screen=True, window_size=(1200, 700))
    p.set_background("white")
    p.add_mesh(ground, color=(0.86, 0.77, 0.58), smooth_shading=True, specular=0.1)
    p.add_mesh(water, scalars="speed", cmap=FOAM, clim=(0, 0.6), opacity=0.9, smooth_shading=True,
               specular=0.4, scalar_bar_args=dict(title="speed (m/s)", color="#374151", n_labels=4,
                                                  position_x=0.78, position_y=0.08, width=0.18, height=0.06))
    p.add_text(f"t = {t:5.1f} s", position="upper_left", font_size=11, color="#374151")
    p.camera_position = [(1.5, -4.5, 9.5), (13.5, 14.2, -1.2), (0, 0, 1)]
    frames.append(Image.fromarray(p.screenshot(return_img=True)))
    p.close()

os.makedirs(os.path.dirname(gif) or ".", exist_ok=True)
small = [f.resize((720, 420), Image.LANCZOS).quantize(colors=128) for f in frames]
small[0].save(gif, save_all=True, append_images=small[1:], duration=180, loop=0, optimize=True)
frames[len(frames) // 2].save(os.path.splitext(gif)[0] + ".png")
print(f"{len(frames)} frames -> {gif}")
