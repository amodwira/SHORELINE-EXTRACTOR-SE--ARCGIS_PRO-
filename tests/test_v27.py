import os
import sys, types, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
arcpy.sa = sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
from skimage import measure

# Scene: Gulf rows 0-14; marsh rows 15-59 at MHW+0.05 (uncertain); a real tidal pond; a sparse mudflat platform.
H=W=60; mhw=0.09
dem = np.full((H,W), 0.20); flat=np.zeros((H,W))
dem[:15,:] = 0.0; flat[:15,:] = 1.0                  # Gulf, hydro-flattened near shore
ndvi = np.full((H,W), 0.55); ndvi[:15,:] = -0.15     # marsh green, water negative
pond = np.zeros((H,W),bool); pond[35:42,10:18]=True
dem[pond]=0.0; flat[pond]=1.0; ndvi[pond]=-0.1
mud = np.zeros((H,W),bool); mud[20:24,40:55]=True; ndvi[mud]=0.12; dem[mud]=0.15   # sparse mud platform >= MHW

# MNDWI result: Gulf water; pond water; BUT wet marsh misread as water in two patches:
spec = np.zeros((H,W),bool); spec[:15,:]=True; spec[pond]=True
wet_edge = np.zeros((H,W),bool); wet_edge[15:19, 0:30]=True     # wet marsh fringe attached to Gulf
wet_int  = np.zeros((H,W),bool); wet_int[45:52, 30:45]=True     # flooded marsh interior
spec |= wet_edge | wet_int | mud

land, water, zone, veg, z = m._dem_zones(dem, flat, mhw, 0.5, 0.5, True, ndvi, 0.3, 0.15)
mv, mg = m._marsh_to_land(zone, ndvi, 0.3, True, z, flat, mhw, use_ground_rule=False)
fused = m._apply_zones(spec, land | mv | mg, water)
assert not fused[wet_edge].any(), "wet marsh fringe -> land (shoreline moves back to Gulf edge)"
assert not fused[wet_int].any(), "flooded marsh interior -> land (no false interior shoreline)"
assert fused[:15,:].all() and fused[pond].all(), "real open water kept"
assert fused[mud].all(), "sparse mudflat stays MNDWI water when ground rule off"
# ground rule on
mv2, mg2 = m._marsh_to_land(zone, ndvi, 0.3, True, z, flat, mhw, use_ground_rule=True)
fused2 = m._apply_zones(spec, land | mv2 | mg2, water)
assert not fused2[mud].any(), "ground>=MHW rule turns mud platform to land"
assert fused2[pond].all(), "flattened pond NOT turned to land by ground rule"

def edge_len(b):
    return sum(len(c) for c in measure.find_contours(b.astype(float), 0.5))
valid=np.ones((H,W),bool)
o_spec,_ = m._ocean_connected(spec, valid)
o_fused,_ = m._ocean_connected(fused, valid)
print("Gulf shoreline row (MNDWI only):", np.where(o_spec[:,5])[0].max(), "| v2.7:", np.where(o_fused[:,5])[0].max(), "| true: 14")
print("water-to-land by veg rule:", int((spec & mv).sum()), "px; zone codes present:", sorted(set(np.unique(zone))))
print("ALL v2.7 TESTS PASSED")
