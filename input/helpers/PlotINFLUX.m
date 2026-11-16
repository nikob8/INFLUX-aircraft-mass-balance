S = shaperead('Indianapolis_Interstate_Mileposts.shp');
mstruct = defaultm('tranmerc');
mstruct.geoid = almanac('earth','grs80','survey foot');
mstruct.origin = [37.5 -85.66666666666667 0];
mstruct.falseeasting = 328083.3333333333;
mstruct.falsenorthing = 820208.3333333331;
mstruct.scalefactor = 0.9999666666666667;
mstruct = defaultm(mstruct);
for i=1:1469; XX(i) = S(i).X; YY(i) = S(i).Y; end
[LAT, LON] = minvtran(mstruct, XX, YY);
scatterm(LAT,LON,30,'.','MarkerEdgeColor',[0.1 0.1 0.1])