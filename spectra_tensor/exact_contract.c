/* Exact Gaussian-integer tensor contractions, accelerated with GMP. */
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
_Static_assert(sizeof(unsigned long)>=sizeof(int64_t), "64-bit coefficient path required");
typedef struct { size_t rows,cols; mpz_t *re,*im; } matrix;
void matrix_free(matrix*m){
    if(!m)return;
    for(size_t i=0;i<m->rows*m->cols;i++){mpz_clear(m->re[i]);mpz_clear(m->im[i]);}
    free(m->re);free(m->im);free(m);
}
matrix* matrix_new(size_t r,size_t c){
    if(!r||!c||r>4096||c>4096||r>SIZE_MAX/c)return NULL;
    matrix*m=calloc(1,sizeof(matrix));if(!m)return NULL;
    m->rows=r;m->cols=c;m->re=malloc(r*c*sizeof(mpz_t));m->im=malloc(r*c*sizeof(mpz_t));
    if(!m->re||!m->im){free(m->re);free(m->im);free(m);return NULL;}
    for(size_t i=0;i<r*c;i++){mpz_init(m->re[i]);mpz_init(m->im[i]);}return m;
}
matrix* matrix_one(void){matrix*m=matrix_new(1,1);if(m)mpz_set_ui(m->re[0],1);return m;}
static void am(mpz_t out,const mpz_t x,int64_t a){
    if(a>0)mpz_addmul_ui(out,x,(unsigned long)a);
    else if(a<0)mpz_submul_ui(out,x,(unsigned long)(-a));
}
static void cm(mpz_t r,mpz_t i,const mpz_t x,const mpz_t y,int64_t a,int64_t b){am(r,x,a);am(r,y,-b);am(i,x,b);am(i,y,a);}
int matrix_add(matrix*a,const matrix*b){
    if(!a||!b||a->rows!=b->rows||a->cols!=b->cols)return 0;
    for(size_t i=0;i<a->rows*a->cols;i++){mpz_add(a->re[i],a->re[i],b->re[i]);mpz_add(a->im[i],a->im[i],b->im[i]);}return 1;
}
matrix* matrix_transfer(const matrix*e,const int64_t*ar,const int64_t*ai,size_t ra,const int64_t*br,const int64_t*bi,size_t rb,int conjugate){
    if(!e)return NULL;
    size_t m=e->rows,n=e->cols;matrix*t=matrix_new(m,rb);if(!t)return NULL;
    for(size_t k=0;k<n;k++)for(size_t j=0;j<rb;j++){
        int64_t r=br[k*rb+j],im=bi[k*rb+j];if(!r&&!im)continue;
        for(size_t i=0;i<m;i++)cm(t->re[i*rb+j],t->im[i*rb+j],e->re[i*n+k],e->im[i*n+k],r,im);
    }
    matrix*out=matrix_new(ra,rb);if(!out){matrix_free(t);return NULL;}
    for(size_t k=0;k<m;k++)for(size_t i=0;i<ra;i++){
        int64_t r=ar[k*ra+i],im=ai[k*ra+i];if(conjugate)im=-im;if(!r&&!im)continue;
        for(size_t j=0;j<rb;j++)cm(out->re[i*rb+j],out->im[i*rb+j],t->re[k*rb+j],t->im[k*rb+j],r,im);
    }
    matrix_free(t);return out;
}
char* matrix_scalar(const matrix*m,int imaginary){
    if(!m||m->rows!=1||m->cols!=1)return NULL;
    const __mpz_struct*x=imaginary?m->im[0]:m->re[0];char*s=malloc(mpz_sizeinbase(x,10)+3);if(s)mpz_get_str(s,10,x);return s;
}
void string_free(char*s){free(s);}
