"""ParaView rendering of the island run: terrain and water surface in relief,
water coloured by speed. Writes one png per output time.

    pvbatch scripts/paraview_island.py results/island_c/island.pvd figures/frames
    python3 scripts/make_gif.py figures/frames figures/island_c.gif

Needs ParaView 5.10 or later (Threshold API). Heights are exaggerated 10x
"""
import os
import sys

from paraview.simple import (Calculator, CellDatatoPointData, ColorBy, CreateView, GetAnimationScene,
                             GetColorTransferFunction, GetScalarBar, PVDReader, Render, SaveScreenshot, Show,
                             Threshold, WarpByScalar)

pvd = sys.argv[1] if len(sys.argv) > 1 else "results/island_c/island.pvd"
out = sys.argv[2] if len(sys.argv) > 2 else "figures/frames"
S = 10
os.makedirs(out, exist_ok=True)

src = PVDReader(FileName=pvd)
pts = CellDatatoPointData(Input=src)

ground = WarpByScalar(Input=pts)
ground.Scalars = ["POINTS", "b"]
ground.ScaleFactor = S

wet = Threshold(Input=pts)
wet.Scalars = ["CELLS", "h"]
wet.LowerThreshold = 1e-3
wet.UpperThreshold = 1e3
water = WarpByScalar(Input=wet)
water.Scalars = ["POINTS", "eta"]
water.ScaleFactor = S
speed = Calculator(Input=water)
speed.ResultArrayName = "speed"
speed.Function = "sqrt(u^2+v^2)"

view = CreateView("RenderView")
view.ViewSize = [1600, 900]
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass  # older ParaView: the background is used as set

g = Show(ground, view)
g.ColorArrayName = ["POINTS", ""]
g.DiffuseColor = [0.86, 0.77, 0.58]
g.AmbientColor = [0.86, 0.77, 0.58]

w = Show(speed, view)
ColorBy(w, ("POINTS", "speed"))
# still water blue, fast water close to white, like foam
lut = GetColorTransferFunction("speed")
lut.RGBPoints = [0, 0.16, 0.47, 0.84, 0.6, 0.92, 0.96, 1]
lut.AutomaticRescaleRangeMode = "Never"
w.Opacity = 0.9
w.Specular = 0.4
w.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = "speed (m/s)"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0.22, 0.25, 0.32]
bar.RangeLabelFormat = "%-#.1f"
bar.WindowLocation = "Upper Right Corner"

# the first render resets the camera, so place it afterwards
Render(view)

view.CameraPosition = [1.5, -4.5, 9.5]
view.CameraFocalPoint = [13.5, 14.2, -1.2]
view.CameraViewUp = [0, 0.0, 1]
view.CameraViewAngle = 30

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
for k, t in enumerate(src.TimestepValues):
    view.ViewTime = t
    Render(view)
    SaveScreenshot(os.path.join(out, f"frame_{k:03d}.png"), view, ImageResolution=[1600, 900])
print(f"{len(src.TimestepValues)} frames in {out}")
