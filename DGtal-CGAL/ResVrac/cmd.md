
# basic:
./tv-triangulation-color -i ../../../depixelize/images/smw_yoshi_input.ppm -b 20 --discontinuities  0.5  -l 0.5 -q 20 --amplitude 0.5 -t 1  -e test.eps -C dual.eps  -E dualMesh.eps --numColorExportEPSDual 2 --displayMesh


# new:

/tv-triangulation-color  -b 20 --discontinuities 1 --amplitude 0.9 -i ../pixelart/mario-2.ppm --fixDarkEdges 20 -R 10 -C reconsMarioOpti.eps  && epstopdf reconsMarioOpti.eps
 
 ./tv-triangulation-color  -b 8 --discontinuities 1 --amplitude 0.9 -i ../pixelart/mario-yoshi.ppm --fixDarkEdges 50 -R 10 -C reconsMarioYoshiOpti.eps  && epstopdf reconsMarioYoshiOpti.eps