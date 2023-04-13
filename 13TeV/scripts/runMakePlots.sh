g++ -Wno-deprecated makePlots.C -o makePlots.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs` -lHist -lCore -lMathCore

echo "************** compiled makePlots.C --- >  makePlots.exe *******************"
echo ">>>>>>>>>>>>>> now running --- makePlots.exe <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
echo ""
./makePlots.exe
echo ""
echo "_______________________ done ______________________________"
rm makePlots.exe
