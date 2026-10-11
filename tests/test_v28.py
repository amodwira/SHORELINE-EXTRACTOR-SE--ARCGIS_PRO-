import os
import sys, types, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
arcpy.sa = sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
from rasterio import features
from rasterio.transform import from_origin
from shapely.geometry import shape, LineString, box
from shapely.ops import unary_union

cs=10.0; H=W=60; T=from_origin(0, H*cs, cs, cs)
ndvi=np.full((H,W),0.6); ndvi[:15,:]=-0.15
ndvi[15:17, 0:30]=0.08     # beach 2 px wide (cols 0-29)
ndvi[15, 30:40]=0.08       # beach 1 px wide (cols 30-39)
# cols 40-59: marsh directly against the Gulf
zone=np.zeros((H,W),np.uint8); zone[:15,:]=2
water=np.zeros((H,W),bool); water[:15,:]=True

marsh,_ = m._marsh_to_land(zone, ndvi, 0.3, True)
mpoly = unary_union([shape(g) for g,v in features.shapes(marsh.astype(np.uint8), transform=T) if v==1])

wpoly = unary_union([shape(g) for g,v in features.shapes(water.astype(np.uint8), transform=T) if v==1])
frame = box(0,0,W*cs,H*cs).exterior
def segments(poly):
    c=list(poly.exterior.coords); segs=[]
    for a,b in zip(c[:-1],c[1:]):
        s=LineString([a,b])
        # split long straight edges into pixel-length pieces (SplitLine works on vertices; densify like a stair-step)
        n=max(1,int(round(s.length/cs)))
        for i in range(n):
            segs.append(LineString([s.interpolate(i/n,True), s.interpolate((i+1)/n,True)]))
    return [s for s in segs if s.distance(frame) > 1e-6 or not frame.buffer(1e-6).contains(s)]
for label, wp in (("unsimplified", wpoly), ("simplified(0.5px)", wpoly.simplify(0.5*cs))):
    segs=[s for s in segments(wp) if not frame.buffer(0.01).contains(s)]
    buf=mpoly.buffer(0.75*cs)
    kept=[s for s in segs if not buf.contains(s.interpolate(0.5,True))]
    excl=[s for s in segs if buf.contains(s.interpolate(0.5,True))]
    xs_k=sorted(set(round(s.interpolate(0.5,True).x/cs-0.5) for s in kept))
    xs_e=sorted(set(round(s.interpolate(0.5,True).x/cs-0.5) for s in excl))
    print(f"{label}: kept cols {xs_k[0]}-{xs_k[-1]} ({len(kept)} segs), excluded cols {xs_e[0]}-{xs_e[-1]} ({len(excl)} segs)")
    assert all(x<40 for x in xs_k) and all(x>=39 for x in xs_e)
print("ALL v2.8 GEOMETRY TESTS PASSED")
