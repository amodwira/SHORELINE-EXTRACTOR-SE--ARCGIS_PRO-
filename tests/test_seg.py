import os
import sys, types, numpy as np
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
arcpy.sa = sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
H=W=20; valid=np.ones((H,W),bool); valid[:, -1]=False
water=np.zeros((H,W),bool); water[:5,:]=True; water[12:15,3:6]=True   # Gulf + inland pond
marsh=np.zeros((H,W),bool); marsh[5:9,10:]=True
ocean,_=m._ocean_connected(water, valid)
seg=m._segmentation(valid, water, ocean, marsh)
assert seg[0,0]==1 and seg[13,4]==2 and seg[18,0]==3 and seg[6,12]==4 and seg[0,-1]==255
print({m.SEG_CLASSES[k][0]: int((seg==k).sum()) for k in (1,2,3,4)})
print("SEGMENTATION TEST PASSED")
