#!/bin/csh

bruk2pipe -in $argv[1] \
  -bad 0.0 -ext -aswap -AMX -decim 2000 -dspfvs 20 -grpdly 67.9862518310547  \
  -xN              2048  -yN                70  \
  -xT              1024  -yT                35  \
  -xMODE            DQD  -yMODE  Echo-AntiEcho  \
  -xSW        10000.000  -ySW         1785.714  \
  -xOBS         500.131  -yOBS         125.763  \
  -xCAR           2.500  -yCAR          43.667  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D         Complex  \
| nmrPipe  -fn SP -off 0.5 -end 1.0 -pow 2 -c 1.0	\
# | nmrPipe  -fn GM -g1 8.0 -g2 18.8			\
| nmrPipe  -fn ZF -size 2048				\
| nmrPipe  -fn FT					\
| nmrPipe  -fn PS -p0 99.5 -p1 0.0 -di			\
# | nmrPipe  -fn POLY -ord 1				\
| nmrPipe  -fn EXT -x1 1.0ppm -xn 6.0ppm  -sw		\
| nmrPipe  -fn TP					\
# | nmrPipe  -fn LP					\
| nmrPipe  -fn SP -off 0.5 -end 1.0 -pow 2 -c 1.0	\
# | nmrPipe  -fn GM -g1 5.0 -g2 13.9			\
| nmrPipe  -fn ZF -size 512				\
# | nmrPipe  -fn FT -neg				\
| nmrPipe  -fn FT					\
| nmrPipe  -fn PS -p0 0.0 -p1 180.0 -di			\
#| nmrPipe  -fn POLY -auto				\
#| nmrPipe  -fn TP					\
#| nmrPipe  -fn POLY -auto				\
#| nmrPipe  -fn EXT -x1 5.5ppm -xn 8.0ppm  -sw		\
#| nmrPipe  -fn TP					\
 -verb -ov -out $argv[2]

# *0.101329118
# e(3.141592653*8*0.08519679836422147072-(0.6*3.141592653*18.8*0.08519679836422147072)^2) < 0.001
# LP: e(3.141592653*5*0.056299997916900077*2-(0.6*3.141592653*13.9*0.056299997916900077*2)^2) < 0.001
#  -xCAR           4.725  -yCAR         122.679  \
#  -xCAR           4.725  -yCAR         132.676  \
