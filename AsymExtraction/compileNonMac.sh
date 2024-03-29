#!/bin/bash
rm PlotRes_Dict.*
CFLAGS="-std=gnu++0x -Wall -ggdb `root-config --cflags --libs` "
#-DDEFAULT_SYSROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -DC_INCLUDE_DIRS=:/usr/local/include:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/11.0.0/include:/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include:/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include:/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks
#CFLAGS=" -Wall -ggdb `root-config --cflags --libs` "
echo cflags: $CFLAGS

rootcint -f PlotRes_Dict.C -c  PlotResults.h

c++ $CFLAGS TwoHadAsymsCleaned.cc TwoHadAsymsCommons.cc MultiPlotter.cxx HadronPairArray.cxx PlotRes_Dict.C PlotResults.cxx -o TwoHadAsymsCMod
#clang $CFLAGS TwoHadAsymsCleaned.cc TwoHadAsymsCommons.cc MultiPlotter.cxx HadronPairArray.cxx PlotRes_Dict.C PlotResults.cxx -o TwoHadAsymsCMod 
c++ $CFLAGS CombPlots.cc TwoHadAsymsCommons.cc PlotRes_Dict.C MultiPlotter.cxx HadronPairArray.cxx PlotResults.cxx -o CombPlots
#c++ $CFLAGS treeTester.cc -o treeTester
#clang $CFLAGS CombPlots.cc TwoHadAsymsCommons.cc PlotRes_Dict.C MultiPlotter.cxx HadronPairArray.cxx PlotResults.cxx -o CombPlots
#c++ $CFLAGS SysRatio.cc TwoHadAsymsCommons.cc PlotRes_Dict.C MultiPlotter.cxx HadronPairArray.cxx PlotResults.cxx -o SysRatio
#c++ $CFLAGS TwoHadAsymsCleaned.cc TwoHadAsymsCommons.cc MultiFitter.cxx HadronQuadArray.cxx HadronPairArrays.cxx -o TwoHadAsymsC
