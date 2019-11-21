#ifndef _GENO_H
#define _GENO_H
#ifdef  __cplusplus
extern "C" {
#endif
	void get_Genotypes(char* chr, int variant_id,int* samples,int numberOfSamples, char* genoFilter);
#ifdef  __cplusplus
}
#endif
#endif
