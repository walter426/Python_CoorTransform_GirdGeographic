from CoorTransform_GirdGeographic import * 

CT_GG = CoorTransformer_GridAndGeographic()

LngLat = CT_GG.CoorTransform_Hk1980ToWgs84(814561, 827359)
EastNorth = CT_GG.CoorTransform_Wgs84ToHK1980(LngLat[0], LngLat[1])

print (LngLat)
print (EastNorth)