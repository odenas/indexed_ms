#ifndef dbwt_h
#define dbwt_h
#ifndef uchar
typedef unsigned char uchar;
#endif
int dbwt_sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs);
int dbwt_sais_int(const int *T, int *SA, int n, int k);
int dbwt_sais(const unsigned char *T, int *SA, int n);
uchar * dbwt_bwt(uchar * T,long n,unsigned int *_last,unsigned int free_text);

#endif
