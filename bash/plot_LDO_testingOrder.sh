gnuplot -e "outputFile='plots/VCC_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VCC';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VCC_LDO_script.gpi
gnuplot -e "outputFile='plots/VEE1_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VEE1';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VEE1_LDO_script.gpi
gnuplot -e "outputFile='plots/VEE_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VEE';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VEE_LDO_script.gpi
gnuplot -e "outputFile='plots/VCC1_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VCC1';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VCC1_LDO_script.gpi
gnuplot -e "outputFile='plots/VFET_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VFET';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VFET_LDO_script.gpi
gnuplot -e "outputFile='plots/VB1_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB1';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VB1_LDO_script.gpi
gnuplot -e "outputFile='plots/VB2_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB2';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VB2_LDO_script.gpi
gnuplot -e "outputFile='plots/VB_testingOrder.eps';\
    xTitleI='Testing order';\
    yTitleI='Current [mA]';\
    xTitleV='Testing order';\
    yTitleV='Voltage [V]';\
    plotTitle='HE board CC4 bias V and I: VB';\

     inputFile='data/Raw_LDO_testingOrder.tsv';" gnuplotScripts/VB_LDO_script.gpi
