gnuplot -e "inputFile='/Users/seacow/Dropbox/seacow_13Inch2014_LBNL/workspace/CMB_S4/feb2024/Run66/Timing/time_14.txt'; ; xTitle='Time []'; yTitle='Voltage [arb]'; outputFile='rawTiming.eps'; plotTitle='Raw Timing Data'; xTitleV='Time [s]'; xTitleI='Time [s]'; yTitleV='Voltage [V]'; yTitleI='SQUID controller Voltage [V]'" timing_raw.gpi
epstopdf rawTiming.eps
