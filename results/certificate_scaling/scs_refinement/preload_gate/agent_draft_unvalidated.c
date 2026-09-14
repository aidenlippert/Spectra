#include <stdio.h>
#include <stdlib.h>
typedef int integer; typedef double doublereal;
static long calls=0,queries=0;
extern void scipy_dsyevd_(char*,char*,integer*,doublereal*,integer*,doublereal*,doublereal*,integer*,integer*,integer*,integer*);
void dsyev_(char *jobz,char *uplo,integer *n,doublereal *a,integer *lda,doublereal *w,doublereal *work,integer *lwork,integer *info){
 calls++; integer li=-1,il=-1,iwq[1]; doublereal wq[1];
 if(*lwork==-1){queries++;scipy_dsyevd_(jobz,uplo,n,a,lda,w,wq,&li,iwq,&il,info);work[0]=wq[0];*info=0;return;}
 scipy_dsyevd_(jobz,uplo,n,a,lda,w,wq,&li,iwq,&il,info); li=(integer)wq[0];il=iwq[0];
 if(li<1||il<1){*info=-8;return;} doublereal *ww=malloc(sizeof(doublereal)*li);integer *iw=malloc(sizeof(integer)*il);
 if(!ww||!iw){free(ww);free(iw);*info=-100;return;} scipy_dsyevd_(jobz,uplo,n,a,lda,w,ww,&li,iw,&il,info);free(ww);free(iw);
}
__attribute__((destructor)) static void done(void){fprintf(stderr,"spectra_evd_shim calls=%ld queries=%ld\n",calls,queries);}
