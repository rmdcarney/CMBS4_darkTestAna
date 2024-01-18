gnuplot -e "outputFile='plots/VCC.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VCC';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VCC_LDO_script.gpi
gnuplot -e "outputFile='plots/VEE1.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VEE1';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VEE1_LDO_script.gpi
gnuplot -e "outputFile='plots/VEE.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VEE';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VEE_LDO_script.gpi
gnuplot -e "outputFile='plots/VCC1.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VCC1';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VCC1_LDO_script.gpi
gnuplot -e "outputFile='plots/VFET.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VFET';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VFET_LDO_script.gpi
gnuplot -e "outputFile='plots/VB1.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB1';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VB1_LDO_script.gpi
gnuplot -e "outputFile='plots/VB2.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB2';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VB2_LDO_script.gpi
gnuplot -e "outputFile='plots/VB.eps';\
    xTitleI='Board #';\
    yTitleI='Current [mA]';\
    xTitleV='Board #';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB';\

     inputFile='data/Raw_LDO.tsv';" gnuplotScripts/VB_LDO_script.gpi
