"""Frames (frame_000.png, ...) to a looping GIF small enough for a README.

    python3 scripts/make_gif.py figures/frames figures/island_c.gif
"""
import glob
import os
import sys

from PIL import Image

src, gif = sys.argv[1], sys.argv[2]
frames = [Image.open(p).convert("RGB") for p in sorted(glob.glob(os.path.join(src, "frame_*.png")))]
w = 760
small = [f.resize((w, round(f.height * w / f.width)), Image.LANCZOS).quantize(colors=128) for f in frames]
small[0].save(gif, save_all=True, append_images=small[1:], duration=180, loop=0, optimize=True)
print(f"{len(frames)} frames -> {gif}, {os.path.getsize(gif) / 1e6:.1f} MB")
