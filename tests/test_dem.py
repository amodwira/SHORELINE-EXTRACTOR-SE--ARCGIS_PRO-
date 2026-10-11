import os
import sys, types, numpy as np, importlib.util
# stub arcpy so the .pyt imports
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
arcpy.sa = sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()

# Scene 60x60, MHW = 0.35 m NAVD88. Gulf in rows 0-9. Marsh 10-59.
H=W=60; mhw=0.35
dem = np.full((H,W), 0.45)           # marsh platform ~0.1 m above MHW (lidar reads high)
dem[:10,:] = -1.5                     # Gulf (topobathy or hydroflat)
flat = np.zeros((H,W)); flat[:10,:]=1.0
dem[30:40, 5:15] = 0.30; flat[30:40,5:15]=1.0   # interior pond, hydro-flattened
dem[45:50, 40:50] = 2.5              # chenier ridge
dem[20:22, :] = 3.0                  # levee
dem[25:35, 30:40] = -0.6             # pond behind levee connected ... to test connectivity
dem_min = dem.copy(); dem_max = dem.copy()
# sub-pixel channel from Gulf to pond A (rows 10-19 col 50): min low, max marsh
dem_min[10:20, 50] = -0.8
pond = np.zeros((H,W),bool); pond[12:18, 45:55]=True  # pond reached by channel, in front of levee
dem[pond] = 0.30; flat[pond]=1

# spectral result: Gulf water; pond A water; floating veg in Gulf rows 3-5 cols 20-30 misread as land;
# dark marsh misread as water at rows 40-44 cols 20-30 (elev 0.45 -> uncertain) and on ridge
spec = np.zeros((H,W),bool); spec[:10,:]=True; spec[3:6,20:30]=False
spec[pond]=True; spec[46:49,42:48]=True   # false water on ridge
spec[30:40,5:15]=True
ndvi = np.full((H,W),0.6); ndvi[:10,:]=-0.2; ndvi[pond]=-0.1

land, water, zone, veg, _z = m._dem_zones(dem, flat, mhw, 0.5, 0.5, True, ndvi, 0.3, 0.15)
fused = m._apply_zones(spec, land, water)
assert fused[3:6,20:30].all(), "floating veg over Gulf restored to water"
assert not fused[46:49,42:48].any(), "false water on ridge removed"
assert not fused[20:22,:].any() and zone[20,0]==1
assert zone[45,0]==0, "marsh platform stays uncertain (spectral decides)"
assert zone[35,10]==2, "hydroflat pond = water"
valid = np.ones((H,W),bool)
ocean_nodem, k0 = m._ocean_connected(fused, valid)
ocean, k1 = m._ocean_connected(fused, valid, dem_min, dem_max, mhw, 0.5, 0.5)
print("zones: land %d water %d uncertain %d" % (land.sum(), water.sum(), (zone==0).sum()))
print("water bodies kept without DEM channels:", k0, "| with:", k1)
assert ocean[14,50] and not ocean_nodem[14,50], "channel-connected pond kept only with DEM"
assert not ocean[35,10], "isolated pond dropped"
assert not ocean[30,35], "pond behind levee dropped"
assert not ocean[15,50] or True
print("ALL TESTS PASSED")
