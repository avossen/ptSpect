#!/bin/bash
rm PlotRes_Dict.*
CFLAGS="-std=gnu++0x -Wall -ggdb `root-config --cflags --libs` "
echo cflags: $CFLAGS
rootcint -f PlotRes_Dict.C -c  PlotResults.h 
c++ $CFLAGS TwoHadAsymsCleaned.cc TwoHadAsymsCommons.cc MultiPlotter.cxx HadronPairArray.cxx PlotRes_Dict.C PlotResults.cxx -o TwoHadAsymsCMod 
c++ $CFLAGS CombPlots.cc TwoHadAsymsCommons.cc PlotRes_Dict.C MultiPlotter.cxx HadronPairArray.cxx PlotResults.cxx -o CombPlots
c++ $CFLAGS SysRatio.cc TwoHadAsymsCommons.cc PlotRes_Dict.C MultiPlotter.cxx HadronPairArray.cxx PlotResults.cxx -o SysRatio
#c++ $CFLAGS TwoHadAsymsCleaned.cc TwoHadAsymsCommons.cc MultiFitter.cxx HadronQuadArray.cxx HadronPairArrays.cxx -o TwoHadAsymsC