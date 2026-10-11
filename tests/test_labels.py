import os
import sys, types, os, tempfile
arcpy = types.ModuleType("arcpy"); sa = types.ModuleType("arcpy.sa"); sa.__all__=[]
calls=[]
class Mgmt:
    def AddColormap(self,p,t,c): calls.append(("clr",open(c).read()))
    def BuildRasterAttributeTable(self,p,o): calls.append(("rat",p))
    def AddField(self,*a,**k): calls.append(("field",a[1]))
rows=[[0,""],[1,""],[2,""],[3,""],[4,""]]
class UC:
    def __init__(s,*a): pass
    def __enter__(s): return iter(rows)
    def __exit__(s,*a): pass
class F:
    def __init__(s,n): s.name=n
arcpy.management=Mgmt(); arcpy.da=types.SimpleNamespace(UpdateCursor=UC); arcpy.ListFields=lambda p:[F("Value")]
arcpy.sa=sa; sys.modules["arcpy"]=arcpy; sys.modules["arcpy.sa"]=sa
from importlib.machinery import SourceFileLoader
m = SourceFileLoader("se", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "Shoreline_Extractor.pyt")).load_module()
# rows updated in-place by cursor? emulate: updateRow is called on cursor; our UC returns iterator without updateRow
class UC2(UC):
    def __enter__(s): return s
    def __iter__(s): return iter(rows)
    def updateRow(s,r): pass
arcpy.da.UpdateCursor=UC2
d=tempfile.mkdtemp(); p=os.path.join(d,"x_dem_zones.tif")
m._label_classified(p, m.ZONE_CLASSES)
print(calls[0][1]); print([r[1] for r in rows])
assert rows[1][1].startswith("High ground") and rows[3][1].startswith("Marsh")
print("LABEL TEST PASSED")
