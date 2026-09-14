#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <dlfcn.h>
typedef int integer; typedef double doublereal;
_Static_assert(sizeof(integer)==4,"LP64 LAPACK requires 32-bit integers");
static unsigned long calls=0,queries=0,failures=0;
typedef void (*evd_fn)(char*,char*,integer*,doublereal*,integer*,doublereal*,doublereal*,integer*,integer*,integer*,integer*);
static evd_fn evd=NULL;
void dsyev_(char *jobz,char *uplo,integer *n,doublereal *a,integer *lda,doublereal *w,doublereal *work,integer *lwork,integer *info){
 *info=0;
 if(*jobz!='V'&&*jobz!='v'&&*jobz!='N'&&*jobz!='n')*info=-1;
 else if(*uplo!='U'&&*uplo!='u'&&*uplo!='L'&&*uplo!='l')*info=-2;
 else if(*n<0)*info=-3;
 else if(*lda<(*n>1?*n:1))*info=-5;
 else if(*lwork!=-1&&(long long)*lwork<(*n>0?3LL*(*n)-1:1))*info=-8;
 if(*info){failures++;return;}
 if(!evd){
   const char *path=getenv("SPECTRA_EVD_LIBRARY");
   void *handle=path?dlopen(path,RTLD_NOW|RTLD_LOCAL|RTLD_DEEPBIND):NULL;
   if(handle)evd=(evd_fn)dlsym(handle,"scipy_dsyevd_");
   if(!evd){*info=-100;failures++;fprintf(stderr,"spectra_evd_shim library load failed\n");return;}
 }
 integer li=-1,il=-1,iwq=0;doublereal wq=0;
 evd(jobz,uplo,n,a,lda,w,&wq,&li,&iwq,&il,info);
 if(*info){failures++;return;}
 if(!isfinite(wq)||wq<1||wq>INT_MAX||iwq<1){*info=-100;failures++;return;}
 work[0]=wq;
 if(*lwork==-1){queries++;return;}
 li=(integer)wq;il=iwq;
 doublereal *ww=malloc(sizeof(doublereal)*(size_t)li);
 integer *iw=malloc(sizeof(integer)*(size_t)il);
 if(!ww||!iw){free(ww);free(iw);*info=-100;failures++;return;}
 calls++;
 evd(jobz,uplo,n,a,lda,w,ww,&li,iw,&il,info);
 if(*info)failures++;
 free(ww);free(iw);
}
__attribute__((destructor)) static void done(void){fprintf(stderr,"spectra_evd_shim real_calls=%lu queries=%lu failures=%lu\n",calls,queries,failures);}
