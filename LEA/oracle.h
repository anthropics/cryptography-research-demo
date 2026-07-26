// Copyright 2026 Anthropic PBC
// SPDX-License-Identifier: Apache-2.0

#ifndef ORACLE_H
#define ORACLE_H
/* Minimal HTTP/1.1 client for the LEA encryption oracle.
   POST /encrypt  {"pt":"<hex>"} -> {"ct":"<hex>"}
   Body cap 1 MiB per message => at most MAX_BLKS blocks per request.
   Oracle address: ORACLE_HOST / ORACLE_PORT environment variables
   (default 127.0.0.1:8666). ORACLE_HOST must be an IPv4 literal.
   The client requires a Content-Length header on responses (no chunked
   transfer encoding) and keep-alive connections, as oracle_server.py
   provides; it exits with an error on any malformed response.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <signal.h>

#ifndef MSG_NOSIGNAL            /* not defined on macOS; SIGPIPE ignored in ora_init */
#define MSG_NOSIGNAL 0
#endif

#define ORACLE_HOST_DEFAULT "127.0.0.1"
#define ORACLE_PORT_DEFAULT 8666
#define MAX_BLKS    32700          /* blocks per request (hex => 32*n bytes body) */

/* ---- optional local simulation mode (for end-to-end testing) ---- */
#include "lea.h"
int ORA_SIM = 0;                 /* if 1, encrypt locally with ORA_SIM_KEY */
uint32_t ORA_SIM_KEY[4];
long long ORA_QUERIES = 0, ORA_BLOCKS = 0;   /* cumulative oracle usage */
static lea_ctx ORA_SIM_CTX;
static int ORA_SIM_READY = 0;

static int ora_fd = -1;

static const char *ora_host(void){ const char*h=getenv("ORACLE_HOST"); return h&&*h? h : ORACLE_HOST_DEFAULT; }
static int ora_port(void){ const char*p=getenv("ORACLE_PORT"); return p&&*p? atoi(p) : ORACLE_PORT_DEFAULT; }

static int ora_connect(void){
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    struct sockaddr_in a; memset(&a,0,sizeof(a));
    a.sin_family=AF_INET; a.sin_port=htons((uint16_t)ora_port());
    if(inet_pton(AF_INET, ora_host(), &a.sin_addr)!=1){
        fprintf(stderr,"ORACLE_HOST must be an IPv4 address: %s\n", ora_host()); exit(4); }
    if(connect(fd,(struct sockaddr*)&a,sizeof(a))<0){ perror("connect"); exit(4); }
    int one=1; setsockopt(fd,IPPROTO_TCP,TCP_NODELAY,&one,sizeof(one));
    return fd;
}
static void ora_init(void){ if(ora_fd<0){ signal(SIGPIPE,SIG_IGN); ora_fd=ora_connect(); } }

static const char HEXD[]="0123456789abcdef";
static inline int hexval(int c){ return (c<='9')? c-'0' : (c|32)-'a'+10; }

static char *REQBUF=NULL, *RESPBUF=NULL;
static size_t REQCAP=0, RESPCAP=0;

/* Encrypt nblk blocks (16 bytes each) from pt -> ct.  nblk <= MAX_BLKS */
static void ora_encrypt(const uint8_t *pt, uint8_t *ct, int nblk){
    ORA_QUERIES++; ORA_BLOCKS += nblk;
    if(ORA_SIM){
        if(!ORA_SIM_READY){ lea_setkey(&ORA_SIM_CTX, ORA_SIM_KEY, 13); ORA_SIM_READY=1; }
        for(int i=0;i<nblk;i++){
            uint32_t x[4];
            for(int k=0;k<4;k++){ const uint8_t*s=pt+16*i+4*k;
                x[k]=(uint32_t)s[0]|((uint32_t)s[1]<<8)|((uint32_t)s[2]<<16)|((uint32_t)s[3]<<24); }
            lea_enc(&ORA_SIM_CTX,x);
            for(int k=0;k<4;k++){ uint8_t*d=ct+16*i+4*k;
                d[0]=x[k];d[1]=x[k]>>8;d[2]=x[k]>>16;d[3]=x[k]>>24; }
        }
        return;
    }
    ora_init();
    size_t hexlen = (size_t)nblk*32;
    size_t need = hexlen + 512;
    if(REQCAP < need){ REQCAP=need*2; REQBUF=realloc(REQBUF,REQCAP);
        if(!REQBUF){fprintf(stderr,"memory allocation failed\n");exit(3);} }
    if(RESPCAP < need){ RESPCAP=need*2; RESPBUF=realloc(RESPBUF,RESPCAP);
        if(!RESPBUF){fprintf(stderr,"memory allocation failed\n");exit(3);} }

    /* build body: {"pt":"<hex>"} */
    char head[256];
    int bodylen = (int)(hexlen + 10);   /* {"pt":"" } => 9 chars + hex ... count precisely below */
    /* body = {"pt":"HEX"}  => 8 + hexlen + 2 = '{','"','p','t','"',':','"',HEX,'"','}' = 7+hexlen+2 = hexlen+9 */
    bodylen = (int)hexlen + 9;
    int hl = snprintf(head,sizeof(head),
        "POST /encrypt HTTP/1.1\r\nHost: %s\r\nContent-Type: application/json\r\nContent-Length: %d\r\nConnection: keep-alive\r\n\r\n",
        ora_host(), bodylen);
    char *p = REQBUF;
    memcpy(p,head,hl); p+=hl;
    memcpy(p,"{\"pt\":\"",7); p+=7;
    for(int i=0;i<nblk*16;i++){ *p++=HEXD[pt[i]>>4]; *p++=HEXD[pt[i]&15]; }
    *p++='"'; *p++='}';
    size_t total = p-REQBUF;

    /* send all (reconnect if the server closed a previous keep-alive) */
    size_t got=0; char *hdr_end=NULL; long clen=-1; size_t hdr_len=0;
    for(int attempt=0; attempt<2; attempt++){
        size_t sent=0; int failed=0;
        while(sent<total){
            ssize_t w = send(ora_fd, REQBUF+sent, total-sent, MSG_NOSIGNAL);
            if(w<=0){ failed=1; break; }
            sent+=w;
        }
        if(!failed){
            got=0; hdr_end=NULL; clen=-1; hdr_len=0;
            while(1){
                if(got+65536 > RESPCAP){ RESPCAP=(got+65536)*2; RESPBUF=realloc(RESPBUF,RESPCAP);
                    if(!RESPBUF){fprintf(stderr,"memory allocation failed\n");exit(3);} }
                ssize_t r = recv(ora_fd, RESPBUF+got, RESPCAP-got-1, 0);
                if(r<=0){ failed=1; break; }
                got+=r; RESPBUF[got]=0;
                if(!hdr_end){
                    hdr_end = strstr(RESPBUF,"\r\n\r\n");
                    if(hdr_end){
                        hdr_len = (hdr_end-RESPBUF)+4;
                        char *cl = strcasestr(RESPBUF,"content-length:");
                        if(cl) clen = strtol(cl+15,NULL,10);
                    }
                }
                if(hdr_end && clen>=0 && got >= hdr_len + (size_t)clen) break;
            }
        }
        if(!failed) break;
        /* reconnect and retry once */
        close(ora_fd); ora_fd = ora_connect();
        if(attempt==1){ fprintf(stderr,"oracle request failed twice\n"); exit(4); }
    }
    /* parse {"ct":"HEX"} */
    char *body = RESPBUF + hdr_len;
    char *q = strstr(body,"\"ct\"");
    if(q){ q = strchr(q+4,':'); if(q) q = strchr(q,'"'); if(q) q++; }
    if(!q || (size_t)(got-(q-RESPBUF)) < (size_t)nblk*32
          || strspn(q,"0123456789abcdefABCDEF") < (size_t)nblk*32){
        fprintf(stderr,"bad oracle response: %.200s\n", body); exit(4); }
    for(int i=0;i<nblk*16;i++){
        ct[i] = (uint8_t)((hexval(q[2*i])<<4) | hexval(q[2*i+1]));
    }
}
#endif
