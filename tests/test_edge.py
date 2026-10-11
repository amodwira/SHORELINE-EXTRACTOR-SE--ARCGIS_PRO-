import os
import sys, types, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
class Point:
    def __init__(s,x,y): s.X=x; s.Y=y
arcpy.Point=Point; arcpy.sa=sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
cs=10.0; H=W=40; grid={"cs":cs,"ll":Point(0.0,0.0)}; ymax=H*cs
valid=np.ones((H,W),bool); valid[-3:,:]=False        # 3 rows of NoData at bottom (inside extent)
# Gulf water rows 25-36; its lower boundary is the NoData edge at row 37 (y = 30)
xs=np.array([(c+0.5)*cs for c in range(W)]); 
edge=m._near_data_edge(xs, np.full(W, ymax-37*cs), valid, grid)      # along NoData edge
beach=m._near_data_edge(xs[5:35], np.full(30, ymax-25*cs), valid, grid)  # real shoreline 12 rows up
assert edge.all() and not beach.any(); print("edge removed:", int(edge.sum()), "| real shoreline kept:", int((~beach).sum()))
print("EDGE TEST PASSED")
