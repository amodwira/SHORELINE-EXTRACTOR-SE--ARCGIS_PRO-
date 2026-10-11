import os
import sys, types, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
class Point:
    def __init__(s,x,y): s.X=x; s.Y=y
arcpy.Point=Point; arcpy.sa=sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
cs=10.0; H=W=60; grid={"cs":cs,"ll":Point(0.0,0.0)}; ymax=H*cs
water=np.zeros((H,W),bool); water[:15,:]=True          # Gulf
pond=np.zeros((H,W),bool); pond[35:42,10:20]=True; water|=pond
marsh=np.zeros((H,W),bool); marsh[15:,:]=True
marsh[15:18, 0:30]=False      # 3-px beach, cols 0-29
marsh[15, 30:40]=False        # 1-px mixed fringe (low NDVI) in front of marsh, cols 30-39
marsh[34,9:21]=marsh[42,9:21]=False; marsh[34:43,9]=marsh[34:43,20]=False   # mixed ring around pond
marsh &= ~water
land=~water
def mid(row_edge, cols):   # midpoints of horizontal edges between row_edge-1 and row_edge
    return np.array([(c+0.5)*cs for c in cols]), np.full(len(cols), ymax-row_edge*cs)
x,y=mid(15, range(60))
old=m._near_marsh(x,y,marsh,grid,0.75); new=m._marsh_majority(x,y,marsh,land,grid,2.5)
print("Gulf edge  old excluded cols:", np.where(old)[0].min(), "| new excluded cols:", np.where(new)[0].min(),"-",np.where(new)[0].max())
assert not new[:30].any(), "3-px beach kept"
assert new[33:].all(), "1-px mixed fringe in front of marsh now excluded (old kept it)"
assert not old[30:39].any()
px,py=mid(42, range(10,20))   # pond southern edge
print("pond edge  old excluded:", int(m._near_marsh(px,py,marsh,grid,0.75).sum()), "/10 | new excluded:", int(m._marsh_majority(px,py,marsh,land,grid,2.5).sum()), "/10")
assert m._marsh_majority(px,py,marsh,land,grid,2.5).all()
print("ALL v2.8.2 TESTS PASSED")
