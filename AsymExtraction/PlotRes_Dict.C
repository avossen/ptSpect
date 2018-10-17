// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PlotRes_Dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "PlotResults.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_PlotResults(void *p = 0);
   static void *newArray_PlotResults(Long_t size, void *p);
   static void delete_PlotResults(void *p);
   static void deleteArray_PlotResults(void *p);
   static void destruct_PlotResults(void *p);
   static void streamer_PlotResults(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PlotResults*)
   {
      ::PlotResults *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PlotResults >(0);
      static ::ROOT::TGenericClassInfo 
         instance("PlotResults", ::PlotResults::Class_Version(), "PlotResults.h", 15,
                  typeid(::PlotResults), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PlotResults::Dictionary, isa_proxy, 16,
                  sizeof(::PlotResults) );
      instance.SetNew(&new_PlotResults);
      instance.SetNewArray(&newArray_PlotResults);
      instance.SetDelete(&delete_PlotResults);
      instance.SetDeleteArray(&deleteArray_PlotResults);
      instance.SetDestructor(&destruct_PlotResults);
      instance.SetStreamerFunc(&streamer_PlotResults);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PlotResults*)
   {
      return GenerateInitInstanceLocal((::PlotResults*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::PlotResults*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PlotResults::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *PlotResults::Class_Name()
{
   return "PlotResults";
}

//______________________________________________________________________________
const char *PlotResults::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotResults*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int PlotResults::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PlotResults*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PlotResults::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotResults*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PlotResults::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PlotResults*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PlotResults::Streamer(TBuffer &R__b)
{
   // Stream an object of class PlotResults.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> meanKinBin1;
      R__b >> meanKinBin2;
      R__b.ReadStaticArray((float*)kTUncertainties);
      R__b.ReadStaticArray((float*)kTUncertainties1);
      R__b.ReadStaticArray((float*)kTUncertainties2);
      R__b.ReadStaticArray((float*)kTSysUncertainties);
      R__b.ReadStaticArray((float*)kTValues);
      R__b.ReadStaticArray((float*)kTValues1);
      R__b.ReadStaticArray((float*)kTValues2);
      R__b.ReadStaticArray((float*)kTMeans);
      R__b >> numKtValues;
      R__b >> exp;
      R__b >> on_res;
      R__b >> isUds;
      R__b >> isCharm;
      R__b >> isMC;
      R__b >> calcType;
      R__b >> binningType;
      R__b >> pidBin;
      R__b >> chargeBin;
      R__b >> firstKinBin;
      R__b >> secondKinBin;
      R__b >> resultIndex;
      R__b >> mPrint;
      R__b.CheckByteCount(R__s, R__c, PlotResults::IsA());
   } else {
      R__c = R__b.WriteVersion(PlotResults::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << meanKinBin1;
      R__b << meanKinBin2;
      R__b.WriteArray(kTUncertainties, 30);
      R__b.WriteArray(kTUncertainties1, 30);
      R__b.WriteArray(kTUncertainties2, 30);
      R__b.WriteArray(kTSysUncertainties, 30);
      R__b.WriteArray(kTValues, 30);
      R__b.WriteArray(kTValues1, 30);
      R__b.WriteArray(kTValues2, 30);
      R__b.WriteArray(kTMeans, 30);
      R__b << numKtValues;
      R__b << exp;
      R__b << on_res;
      R__b << isUds;
      R__b << isCharm;
      R__b << isMC;
      R__b << calcType;
      R__b << binningType;
      R__b << pidBin;
      R__b << chargeBin;
      R__b << firstKinBin;
      R__b << secondKinBin;
      R__b << resultIndex;
      R__b << mPrint;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PlotResults(void *p) {
      return  p ? new(p) ::PlotResults : new ::PlotResults;
   }
   static void *newArray_PlotResults(Long_t nElements, void *p) {
      return p ? new(p) ::PlotResults[nElements] : new ::PlotResults[nElements];
   }
   // Wrapper around operator delete
   static void delete_PlotResults(void *p) {
      delete ((::PlotResults*)p);
   }
   static void deleteArray_PlotResults(void *p) {
      delete [] ((::PlotResults*)p);
   }
   static void destruct_PlotResults(void *p) {
      typedef ::PlotResults current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_PlotResults(TBuffer &buf, void *obj) {
      ((::PlotResults*)obj)->::PlotResults::Streamer(buf);
   }
} // end of namespace ROOT for class ::PlotResults

namespace {
  void TriggerDictionaryInitialization_PlotRes_Dict_Impl() {
    static const char* headers[] = {
"PlotResults.h",
0
    };
    static const char* includePaths[] = {
"/Applications/root_v6.14.04/include",
"/Users/avossen/Documents/myProjects/belle/AsymExtraction/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "PlotRes_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$PlotResults.h")))  PlotResults;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "PlotRes_Dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "PlotResults.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"PlotResults", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("PlotRes_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_PlotRes_Dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_PlotRes_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_PlotRes_Dict() {
  TriggerDictionaryInitialization_PlotRes_Dict_Impl();
}
