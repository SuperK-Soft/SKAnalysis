/********************************************************************
* NTagDataModelDict.h
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************************/
#ifdef __CINT__
#error NTagDataModelDict.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#define G__PRIVATE_GVALUE
#include "G__ci.h"
#include "FastAllocString.h"
extern "C" {
extern void G__cpp_setup_tagtableNTagDataModelDict();
extern void G__cpp_setup_inheritanceNTagDataModelDict();
extern void G__cpp_setup_typetableNTagDataModelDict();
extern void G__cpp_setup_memvarNTagDataModelDict();
extern void G__cpp_setup_globalNTagDataModelDict();
extern void G__cpp_setup_memfuncNTagDataModelDict();
extern void G__cpp_setup_funcNTagDataModelDict();
extern void G__set_cpp_environmentNTagDataModelDict();
}


#include "TObject.h"
#include "TMemberInspector.h"
#include "Candidate.h"
#include "Cluster.h"
#include "EventCandidates.h"
#include "EventParticles.h"
#include "EventTrueCaptures.h"
#include "Particle.h"
#include "PMTHitCluster.h"
#include "PMTHit.h"
#include "TrueCapture.h"
#include <algorithm>
namespace std { }
using namespace std;

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__NTagDataModelDictLN_TClass;
extern G__linked_taginfo G__NTagDataModelDictLN_TBuffer;
extern G__linked_taginfo G__NTagDataModelDictLN_TMemberInspector;
extern G__linked_taginfo G__NTagDataModelDictLN_string;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_maplEstringcOfloatcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOfloatgRsPgRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_Candidate;
extern G__linked_taginfo G__NTagDataModelDictLN_ClusterlECandidategR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlECandidatecOallocatorlECandidategRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlECandidatecOallocatorlECandidategRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlEstringcOallocatorlEstringgRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlEstringcOallocatorlEstringgRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_maplEstringcOvectorlEfloatcOallocatorlEfloatgRsPgRmUcOlesslEstringgRcOallocatorlEpairlEconstsPstringcOvectorlEfloatcOallocatorlEfloatgRsPgRmUgRsPgRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTBaselEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTBaselEdoublegR;
extern G__linked_taginfo G__NTagDataModelDictLN_TVectorTlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TVectorTlEdoublegR;
extern G__linked_taginfo G__NTagDataModelDictLN_TElementActionTlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TElementPosActionTlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTRow_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTRowlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTDiag_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTColumn_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTFlat_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSub_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSparseRow_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSparseDiag_constlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTColumnlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTDiaglEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTFlatlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSublEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSparseRowlEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TMatrixTSparseDiaglEfloatgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TVector3;
extern G__linked_taginfo G__NTagDataModelDictLN_Particle;
extern G__linked_taginfo G__NTagDataModelDictLN_EventParticles;
extern G__linked_taginfo G__NTagDataModelDictLN_ClusterlEParticlegR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlEParticlecOallocatorlEParticlegRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlEParticlecOallocatorlEParticlegRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_TrueCapture;
extern G__linked_taginfo G__NTagDataModelDictLN_EventTrueCaptures;
extern G__linked_taginfo G__NTagDataModelDictLN_ClusterlETrueCapturegR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlETrueCapturecOallocatorlETrueCapturegRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlETrueCapturecOallocatorlETrueCapturegRsPgRcLcLiteratorgR;
extern G__linked_taginfo G__NTagDataModelDictLN_PMTHit;
extern G__linked_taginfo G__NTagDataModelDictLN_ClusterlEPMTHitgR;
extern G__linked_taginfo G__NTagDataModelDictLN_vectorlEPMTHitcOallocatorlEPMTHitgRsPgR;
extern G__linked_taginfo G__NTagDataModelDictLN_reverse_iteratorlEvectorlEPMTHitcOallocatorlEPMTHitgRsPgRcLcLiteratorgR;

/* STUB derived class for protected member access */
class G__NTagDataModelDictLN_EventParticles_PR : public EventParticles {
 public:
};
class G__NTagDataModelDictLN_EventTrueCaptures_PR : public EventTrueCaptures {
 public:
};
typedef Cluster<Candidate> G__ClusterlECandidategR;
typedef Cluster<Particle> G__ClusterlEParticlegR;
typedef Cluster<TrueCapture> G__ClusterlETrueCapturegR;
typedef Cluster<PMTHit> G__ClusterlEPMTHitgR;