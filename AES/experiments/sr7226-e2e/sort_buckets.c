// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <omp.h>
static int cmp_u64(const void *x,const void *y){ uint64_t a=*(const uint64_t*)x,b=*(const uint64_t*)y; return (a>b)-(a<b); }
int main(int argc,char**argv){
    const char *dir=argv[1]; int nthr=argc>2?atoi(argv[2]):32;
    omp_set_num_threads(nthr);
    #pragma omp parallel for schedule(dynamic,1)
    for(int b=0;b<256;b++){
        char in[1024],out[1024]; snprintf(in,sizeof in,"%s/b%02x.bin",dir,b); snprintf(out,sizeof out,"%s/b%02x.sorted",dir,b);
        int fd=open(in,O_RDONLY); if(fd<0){perror(in);exit(1);} off_t sz=lseek(fd,0,SEEK_END); uint64_t n=sz/8;
        uint64_t *a=malloc(sz?sz:8); if(!a){fprintf(stderr,"malloc %lld fail\n",(long long)sz);exit(1);}
        uint64_t off=0; while(off<(uint64_t)sz){ ssize_t r=pread(fd,(char*)a+off,sz-off,off); if(r<=0){perror("pread");exit(1);} off+=r; }
        close(fd);
        qsort(a,n,8,cmp_u64);
        uint64_t m=0; for(uint64_t i=0;i<n;i++) if(i==0||a[i]!=a[m-1]) a[m++]=a[i];   /* dedup */
        FILE*f=fopen(out,"wb"); fwrite(a,8,m,f); fclose(f); free(a); unlink(in);
        fprintf(stderr,"[sort] b%02x: %llu -> %llu unique\n",b,(unsigned long long)n,(unsigned long long)m);
    }
    return 0;
}
