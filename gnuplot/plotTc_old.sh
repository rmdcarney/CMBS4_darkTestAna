gnuplot -e "inputFile='Tc_04_conv.txt'; secondFile='Tc_05_conv.txt'; xTitle='Temperature [K]'; yTitle='Resistance [arb]'; outputFile='Tc_12_13.eps'; plotTitle='R135a-1: Pixel 155 90T, 3mK/min sweep'" Tc.gpi
epstopdf Tc_12_13.eps
#gnuplot -e "inputFile='Tc_12_conv.txt'; secondFile='Tc_13_conv.txt'; xTitle='Temperature [K]'; yTitle='Resistance [arb]'; outputFile='Tc_12_science.eps'; plotTitle='R135a-1: Pixel 155 90T'" Tc_science.gpi
#epstopdf Tc_12_science.eps
