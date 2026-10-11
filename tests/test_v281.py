import os
import sys, types, time, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
class Point:
    def __init__(s,x,y): s.X=x; s.Y=y
arcpy.Point=Point; arcpy.sa = sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
from shapely.geometry import box, Point as SP

# 1) _near_marsh matches exact geometric buffer test on random points
rng=np.random.default_rng(0); H,W,cs=50,60,10.0
grid={"cs":cs,"ll":Point(1000.0,2000.0)}
marsh=rng.random((H,W))<0.08
xs=1000+rng.random(4000)*W*cs; ys=2000+rng.random(4000)*H*cs
# include points exactly on pixel edges (like real segment midpoints)
xs[:1000]=1000+np.round(xs[:1000]/cs)*0+ (rng.integers(0,W,1000))*cs; ys[:1000]=2000+(rng.integers(0,H,1000)+0.5)*cs
fast=m._near_marsh(xs,ys,marsh,grid,0.75)
ymax=2000+H*cs
cells=[box(1000+c*cs, ymax-(r+1)*cs, 1000+(c+1)*cs, ymax-r*cs) for r,c in zip(*np.where(marsh))]
from shapely.strtree import STRtree
tree=STRtree(cells)
exact=np.array([len(tree.query(SP(x,y).buffer(0.75*cs+1e-9, quad_segs=64), predicate="intersects"))>0 for x,y in zip(xs,ys)])
mism=int((fast!=exact).sum()); print("mismatches vs exact geometry:", mism, "of", len(xs))
assert mism <= 3   # only circle-polygonisation tolerance at exact tangency

# 2) speed at real scene size (3871 x 4442) with 300k segments
big=np.random.default_rng(1).random((4442,3871))<0.18
g2={"cs":10.0,"ll":Point(456480.0,3291080.0)}
px=456480+np.random.default_rng(2).random(300000)*38710; py=3291080+np.random.default_rng(3).random(300000)*44420
t=time.time(); f=m._near_marsh(px,py,big,g2,0.75); print(f"300k segments: {time.time()-t:.2f} s, {f.mean()*100:.1f}% near marsh")

# 3) flat vegetated surface is no longer water; flat bare surface still is
dem=np.full((3,3),0.2); flat=np.ones((3,3)); ndvi=np.array([[0.6]*3,[0.05]*3,[np.nan]*3])
land,water,zone,veg,z=m._dem_zones(dem,flat,0.09,0.5,0.5,True,None,0.3,0.0,flat_ndvi=ndvi,flat_ndvi_max=0.3)
assert not water[0].any() and water[1].all() and water[2].all(); print("vegetated-flat guard OK")
print("ALL v2.8.1 TESTS PASSED")
