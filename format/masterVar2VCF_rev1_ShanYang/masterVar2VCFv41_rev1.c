/* This program convert masterVar file to VCF file v4.1.

The input format is masterVar file. It has some lines for meta-informatio and one headline. The data lines follow.
There could be empty fields (two or more "\t" together). 

#ASSEMBLY_ID    GS19240-1100-37-ASM
#.....
#TYPE   VAR-OLPL

>locus	ploidy	chromosome	begin	end	zygosity	varType	reference	allele1Seq	allele2Seq	allele1Score	allele2Score	allele1HapLink	allele2HapLink	xRef	evidenceIntervalId	allele1ReadCount	allele2ReadCount	neitherAlleleReadCount	totalReadCount	allele1Gene	allele2Gene	miRBaseId	repeatMasker	segDupOverlap	relativeCoverage	calledPloidy
21232043	2	chr22	16052238	16052239	het-ref	snp	A	G	A	89	51	3953118	3953117	dbsnp.100:rs2843214;dbsnp.116:rs6518413	13	17	7	0	24				MER4E:ERV1:21.2	1	1.65	N

The output will be in VCF file (v4.1) that contains all the information in the masterVar file.
*/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <time.h>

FILE *srR;
long R_pos = 0;

char get_leading_base(char ref_chr[2000], long pos) {
	char cc; 

	cc = '0';
	while ((R_pos < pos)&&(!feof(srR))) {
	    cc = fgetc (srR);	
	    if ( ((cc>='A')&&(cc<='Z'))||((cc>='a')&&(cc<='z'))) {
		R_pos ++;
		if ((cc>='a')&&(cc<='z')) cc = cc+'A'-'a';
	    }
	
	}
	return (cc);
}

int main(int argc, char **argv) {

        FILE *srV, *des;
	char in_fl[2000], cmd[2000];
        char bufV[50000], s_temp[50000], bufR[50000], pchr[2000];
	char asmID[2000], build[2000], dbsnp[2000];
	time_t  currtime;                                                    
	char charTime[10] = {0};
	int ii, jj, ct, kk, kk0;
	long bg, pbg;
	int a1Sc, a2Sc, QUAL, a1RC, a2RC, tRC;
	char chr[50], zg[50], ref[2000], a1Seq[2000], a2Seq[2000], a1HL[2000], a2HL[2000], xRef[2000];
	char a1G[2000], a2G[2000], mRNA[2000], rmsk[2000], segDup[2000], rCov[2000], cPd[2000];
	int c_chr, c_bg, c_zg, c_ref, c_a1Seq, c_a2Seq, c_a1Sc, c_a2Sc, c_a1HL, c_a2HL, c_xRef, c_a1RC, c_a2RC, c_tRC;
	int c_a1G, c_a2G, c_mRNA, c_rmsk, c_segDup, c_rCov, c_cPd;
	int f_chr, f_bg, f_zg, f_ref, f_a1Seq, f_a2Seq, f_a1Sc, f_a2Sc, f_a1HL, f_a2HL, f_xRef, f_a1RC, f_a2RC, f_tRC;
	int f_a1G, f_a2G, f_mRNA, f_rmsk, f_segDup, f_rCov, f_cPd;
	char FORMAT[5000], F1[200000], Fa[20000];
	char a1Gt, a2Gt, ld_base;
 
        if ( argc != 4 ) {
                printf ("Usage: masterVar2VCFv41 input_masterVar(tsv or tsv.bz2) reference_file(fa.bz2) output_VCF_file\n");
                exit (1);
        }

	strcpy (in_fl, argv[1]);
	if (strstr(in_fl, ".bz2")) {
          snprintf (cmd, 2000, "bzcat %s", argv[1]);
          srV = popen(cmd, "r");
        }
        else srV = fopen(in_fl, "r");
        if (srV == NULL) {
                printf ("%s\n", argv[1]);
                perror("fopen input masterVar file");
                exit(1);
        }

        snprintf (cmd, 2000, "bzcat %s", argv[2]);
        srR = popen(cmd, "r");
        if (srV == NULL) {
                printf ("%s\n", argv[2]);
                perror("fopen input reference file");
                exit(1);
        }

        while ( fgets(bufV, 50000, srV) ) {
	    if (bufV[0]!='#') break;
	    if (strstr(bufV, "ASSEMBLY_ID")) sscanf (bufV, "%*s %s", asmID); 
	    if (strstr(bufV, "GENOME_REFERENCE")) sscanf (bufV, "%*s %*s %*s %s", build); 
	    if (strstr(bufV, "DBSNP_BUILD")) sscanf (bufV, "%*s %*s %*s %s", dbsnp); 
	}
	
        des = fopen(argv[3], "w");
        if (des == NULL) {
                printf ("%s\n", argv[3]);
                perror("fopen output VCF file");
                exit(1);
        }

//// Print meta information
	fprintf (des, "##fileformat=VCFv4.1\n");
	fprintf (des, "##fileDate=");
	time(&currtime);                                                     
	strftime(charTime,sizeof(charTime)-1,"%Y%m%d",localtime(&currtime)); 
	fprintf (des, "%s\n",charTime);  
	fprintf (des, "##source=masterVar2VCFv40\n");
	fprintf (des, "##reference=%s\n", argv[2]);
	for (ii = 1; ii < 26; ii ++ ){
	    if (ii<23) fprintf (des, "##contig=<ID=%d,assembly=B%s,species=\"Homo sapiens\",taxonomy=x>\n", ii, build);
	    if (ii==23) fprintf (des, "##contig=<ID=X,assembly=B%s,species=\"Homo sapiens\",taxonomy=x>\n", build);
	    if (ii==24) fprintf (des, "##contig=<ID=Y,assembly=B%s,species=\"Homo sapiens\",taxonomy=x>\n", build);
	    if (ii==25) fprintf (des, "##contig=<ID=M,assembly=B%s,species=\"Homo sapiens\",taxonomy=x>\n", build);
	}
	fprintf (des, "##phasing=partial\n");

//// Below are the INFO & FORMAT fields
	fprintf (des, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
	fprintf (des, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
//	fprintf (des, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
//	fprintf (des, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n");
	fprintf (des, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build %s\">\n", dbsnp);
//	fprintf (des, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n");
//	fprintf (des, "##FILTER=<ID=q10,Description=\"Quality below 10\">\n");
//	fprintf (des, "##FILTER=<ID=s50,Description=\"Less than 50%% of samples have data\">\n");
	fprintf (des, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf (des, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
	fprintf (des, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
	fprintf (des, "##FORMAT=<ID=HDP,Number=2,Type=Integer,Description=\"Haplotype Read Depth\">\n");
	fprintf (des, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n");
	fprintf (des, "##FORMAT=<ID=PS,Number=2,Type=Integer,Description=\"Phase Set\">\n");
	fprintf (des, "##FORMAT=<ID=GENE,Number=.,Type=String,Description=\"Overlaping Genes\">\n");
	fprintf (des, "##FORMAT=<ID=mRNA,Number=.,Type=String,Description=\"Overlaping mRNA\">\n");
	fprintf (des, "##FORMAT=<ID=rmsk,Number=.,Type=String,Description=\"Overlaping Repeats\">\n");
	fprintf (des, "##FORMAT=<ID=segDup,Number=.,Type=String,Description=\"Overlaping segmentation duplication\">\n");
	fprintf (des, "##FORMAT=<ID=rCov,Number=1,Type=Float,Description=\"relative Coverage\">\n");
	fprintf (des, "##FORMAT=<ID=cPd,Number=1,Type=String,Description=\"called Ploidy(level)\">\n");

//// Print the headline

	fprintf (des, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", asmID);

	a1Sc=a2Sc=QUAL=a1RC=a2RC=tRC=-1;
	c_chr=c_bg=c_zg=c_ref=c_a1Seq=c_a2Seq=c_a1Sc=c_a2Sc=c_a1HL=c_a2HL=c_xRef=c_a1RC=c_a2RC=c_tRC=-1;
	c_a1G=c_a2G=c_mRNA=c_rmsk=c_segDup=c_rCov=c_cPd=-1;
	bg = pbg = -1;

	a1Gt = a2Gt = 'N';

	while (fgets(bufR, 50000, srR)) {
	    if (bufR[0]=='>') break;
	}
	kk = 0;
	while ( (bufR[kk+4]!='\t')&&(bufR[kk+4]!=' ')&&(bufR[kk+4]!='\n')&&(bufR[kk+4]!='\0')) {pchr[kk]=bufR[kk+4]; kk ++;}
	pchr[kk]='\0';

        while ( fgets(bufV, 50000, srV) ) {
	  if (strstr(bufV, "chromosome")) {
	      ii = jj = ct = 0;
	      while (bufV[ii]!='\0') {
		s_temp[jj] = bufV[ii];
		if ((bufV[ii]=='\t')||(bufV[ii]=='\n')) { 
		    s_temp[jj] = '\0'; 
		    jj = -1; 
		    ct ++;
		    if (strcmp(s_temp, "chromosome")==0) {c_chr = ct;}
		    if (strcmp(s_temp, "begin")==0) {c_bg = ct;}
		    if (strcmp(s_temp, "zygosity")==0) {c_zg = ct;}
		    if (strcmp(s_temp, "reference")==0) {c_ref = ct;}
		    if (strcmp(s_temp, "allele1Seq")==0) {c_a1Seq = ct;}
		    if (strcmp(s_temp, "allele2Seq")==0) {c_a2Seq = ct;}
		    if (strcmp(s_temp, "allele1Score")==0) {c_a1Sc = ct;}
		    if (strcmp(s_temp, "allele2Score")==0) {c_a2Sc = ct;}
		    if (strcmp(s_temp, "allele1HapLink")==0) {c_a1HL = ct;}
		    if (strcmp(s_temp, "allele2HapLink")==0) {c_a2HL = ct;}
		    if (strcmp(s_temp, "xRef")==0) {c_xRef = ct;}
		    if (strcmp(s_temp, "allele1ReadCount")==0) {c_a1RC = ct;}
		    if (strcmp(s_temp, "allele2ReadCount")==0) {c_a2RC = ct;}
		    if (strcmp(s_temp, "totalReadCount")==0) {c_tRC = ct;}
		    if (strcmp(s_temp, "allele1Gene")==0) {c_a1G = ct;}
		    if (strcmp(s_temp, "allele2Gene")==0) {c_a2G = ct;}
		    if (strcmp(s_temp, "miRBaseId")==0) {c_mRNA = ct;}
		    if (strcmp(s_temp, "repeatMasker")==0) {c_rmsk = ct;}
		    if (strcmp(s_temp, "segDupOverlap")==0) {c_segDup = ct;}
		    if (strcmp(s_temp, "relativeCoverage")==0) {c_rCov = ct;}
		    if (strcmp(s_temp, "calledPloidy")==0) {c_cPd = ct;}
		}
		ii ++;
		jj ++;
	      }
	      break;
	  }
	}

        while ( fgets(bufV, 50000, srV) ) {
	    ref[0]=a1Seq[0]=a2Seq[0]=a1HL[0]=a2HL[0]=xRef[0]=a1G[0]=a2G[0]=mRNA[0]=rmsk[0]=segDup[0]=rCov[0]=cPd[0]='\0';
	    a1Sc = a2Sc = QUAL = -1;
	    f_chr=f_bg=f_zg=f_ref=f_a1Seq=f_a2Seq=f_a1Sc=f_a2Sc=f_a1HL=f_a2HL=f_xRef=f_a1RC=f_a2RC=f_tRC=0;
	    f_a1G=f_a2G=f_mRNA=f_rmsk=f_segDup=f_rCov=f_cPd=0;

	    ii = jj = ct = 0;
	    while (bufV[ii]!='\0') {
		s_temp[jj] = bufV[ii];
		if ((bufV[ii]=='\t')||(bufV[ii]=='\n')) { 
		    s_temp[jj] = '\0'; 
		    jj = -1; 
		    ct++; 
		    if (ct == c_chr) {
			kk = 0; 
			while (s_temp[kk+3]!='\0') {
			    chr[kk]=s_temp[kk+3]; 
			    kk++;
			} 
			chr[kk]='\0';
		    }	
		    if (c_bg == ct) { bg = atol(s_temp); if (strlen(s_temp)>0) f_bg=1;}
		    if (c_zg == ct) { strcpy(zg, s_temp); if (strlen(s_temp)>0) f_zg=1;}
		    if (c_ref == ct) {strcpy (ref, s_temp); if (strlen(s_temp)>0) f_ref=1;}
		    if (c_a1Seq == ct) {strcpy (a1Seq, s_temp); if (strlen(s_temp)>0) f_a1Seq=1;}
		    if (c_a2Seq == ct) {strcpy (a2Seq, s_temp); if (strlen(s_temp)>0) f_a2Seq=1;}
		    if (c_a1Sc == ct) {a1Sc=atoi(s_temp); if (strlen(s_temp)>0) f_a1Sc=1;}
		    if (c_a2Sc == ct) {a2Sc=atoi(s_temp); if (strlen(s_temp)>0) f_a2Sc=1;}
		    if (c_a1HL == ct) {strcpy (a1HL, s_temp); if (strlen(s_temp)>0) f_a1HL=1;}
		    if (c_a2HL == ct) {strcpy (a2HL, s_temp); if (strlen(s_temp)>0) f_a2HL=1;}
		    if (c_xRef == ct) {strcpy (xRef, s_temp); if (strlen(s_temp)>0) f_xRef=1;}
		    if (c_a1RC == ct) {a1RC=atoi(s_temp); if (strlen(s_temp)>0) f_a1RC=1;}
		    if (c_a2RC == ct) {a2RC=atoi(s_temp); if (strlen(s_temp)>0) f_a2RC=1;}
		    if (c_tRC == ct) {tRC=atoi(s_temp); if (strlen(s_temp)>0) f_tRC=1;}
		    if (c_a1G == ct) {strcpy (a1G, s_temp); if (strlen(s_temp)>0) f_a1G=1;}
		    if (c_a2G == ct) {strcpy (a2G, s_temp); if (strlen(s_temp)>0) f_a2G=1;}
		    if (c_mRNA == ct) {strcpy (mRNA, s_temp); if (strlen(s_temp)>0) f_mRNA=1;}
		    if (c_rmsk == ct) {strcpy (rmsk, s_temp); if (strlen(s_temp)>0) f_rmsk=1;}
		    if (c_segDup == ct) {strcpy (segDup, s_temp); if (strlen(s_temp)>0) f_segDup=1;}
		    if (c_rCov == ct) {strcpy (rCov, s_temp); if (strlen(s_temp)>0) f_rCov=1;}
		    if (c_cPd == ct) {strcpy (cPd, s_temp); if (strlen(s_temp)>0) f_cPd=1;}

		    if (f_a1Sc&&f_a2Sc) {
		    	if (a1Sc < a2Sc) QUAL = a1Sc;
		    	else QUAL = a2Sc; 
		    }
		    else {
			if (f_a1Sc&&(!f_a2Sc)) QUAL = a1Sc;
			if ((!f_a1Sc)&&f_a2Sc) QUAL = a1Sc;
			if ((!f_a1Sc)&&(!f_a2Sc)) QUAL = -1;
		    }

		

		}
		ii ++;
		jj ++;
	    }

	    if ( strcmp(ref, "=") ) {                                //if the entry is not hom ref or complete no-call
		if ((!((strlen(ref)==1)&&(strlen(a1Seq)==1)&&(strlen(a2Seq)==1))&&(strcmp(zg,"hap")))||(!((strlen(ref)==1)&&(strlen(a1Seq)==1))&&(strcmp(zg,"hap")==0))) {	//if need to fetch reference base, not snp
		    
		    while (strcmp(pchr, chr)) {		//new chromosome
			while (fgets(bufR, 50000, srR)) {
			    if (bufR[0]=='>') break;
			}
		    	kk = 0;
			while ( (bufR[kk+4]!='\t')&&(bufR[kk+4]!=' ')&&(bufR[kk+4]!='\n')&&(bufR[kk+4]!='\0')) {pchr[kk]=bufR[kk+4]; kk ++;}
		    	pchr[kk]='\0';
			R_pos = 0;
			pbg = -1;
	
			if (feof(srR)) break;
		    }

		    if ( ((bg>0)&&((ref[0]!=a1Seq[0])||(ref[0]!=a2Seq[0]))&&(strcmp(zg,"hap")))||((bg>0)&&(ref[0]!=a1Seq[0])&&(strcmp(zg,"hap")==0)) ) {//only need to fetch if the leading bases are different
		    	ld_base = get_leading_base (chr, bg);
			if (ld_base !='0') {
		    	  fprintf (des, "%s\t%ld", chr, bg);   
		    	  snprintf (s_temp, 50000, "%c%s", ld_base, ref);
		    	  strcpy (ref, s_temp);
		    	  snprintf (s_temp, 50000, "%c%s", ld_base, a1Seq);
		    	  strcpy (a1Seq, s_temp);
		    	  snprintf (s_temp, 50000, "%c%s", ld_base, a2Seq);
		    	  strcpy (a2Seq, s_temp);
			}
			else fprintf (des, "%s\t%ld", chr, bg+1);
		    }
		    else {
			fprintf (des, "%s\t%ld", chr, bg+1);   
		    }
		}
		else fprintf (des, "%s\t%ld", chr, bg+1);   

		if (f_xRef) {			//has dbsnp annotation
		    kk = kk0 = 0;

		    while (xRef[kk]!='\0') {
			while ((xRef[kk]!=':')&&(xRef[kk]!='\0')) kk ++;
			if (xRef[kk]=='\0') break;
			kk ++;
			while ((xRef[kk]!=';')&&(xRef[kk]!='\0')) {
			    xRef[kk0]=xRef[kk];
			    kk0++;
			    kk++;
			}
			if (xRef[kk]==';') xRef[kk0]=',';
			if (xRef[kk]=='\0') {xRef[kk0]='\0'; break;}
			kk ++;
			kk0 ++;
		    }
		}
		if (f_xRef) fprintf (des, "\t%s", xRef);   
		else fprintf (des, "\t.");   
		
		if (strcmp (zg, "hap")==0 ) {				//haploid 
		  if (strchr(a1Seq,'?')) {				//nocalls in a1 only
			a1Gt='.';
			fprintf (des, "\t%s\t%s", ref, a2Seq); 
			a2Gt='A';
		  }
		  else {
		    if ( strcmp (a1Seq,"N")==0 ) {    //N in the alleles in snp
			fprintf (des, "\t%s\t.", ref); 
			a1Gt='.';
			a2Gt='A';
		    }
		    else {
			fprintf (des, "\t%s\t%s", ref, a1Seq);
			a1Gt = '1';
			a2Gt = 'A';
		    }
		  }
		}

		else {							//not haploid, diploid region
		if (strchr(a1Seq, '?')||strchr(a2Seq, '?')) {		//nocalls
		  if (strchr(a1Seq,'?')&&strchr(a2Seq,'?')) {		//nocalls in both
			a1Gt=a2Gt='.';
			fprintf (des, "\t%s\t.", ref);
		  }
		  else if (strchr(a1Seq,'?')) {				//nocalls in a1 only
			a1Gt='.';
			if (strcmp(a2Seq,ref)==0) { fprintf (des, "\t%s\t.", ref); a2Gt='0';}
			else { fprintf (des, "\t%s\t%s", ref, a2Seq); a2Gt='1';}
		  }
		  else {						//nocalls in a2 only
			a2Gt='.';
			if (strcmp(a1Seq,ref)==0) { fprintf (des, "\t%s\t.", ref); a1Gt='0';}
			else { fprintf (des, "\t%s\t%s", ref, a1Seq); a1Gt='1';}
		  }
		}
		else {						//no ? in the sequence
		 if ( (strcmp (a1Seq,"N")==0)||(strcmp (a2Seq,"N")==0)) {    //N in one of the alleles in snp
		    if (strcmp (a2Seq,"N")==0) {
			fprintf (des, "\t%s\t.", ref); 
			a1Gt='0';
			a2Gt='.';
		    }
		    else {
			fprintf (des, "\t%s\t.", ref); 
			a1Gt='.';
			a2Gt='0';
		    }
		 }
		 else {							//no N or ? in the alleles
		  if (strcmp (a1Seq,a2Seq)==0) {		    //only one ALT genotype
		    if (strcmp (ref, a1Seq)==0) {              //hom ref
			fprintf (des, "\t%s\t.", ref);
			a1Gt = a2Gt = '0';
		    }
		    else {					//hom alt
			fprintf (des, "\t%s\t%s", ref, a1Seq);
			a1Gt = a2Gt = '1';
		    }
		  }
		  else {
		    if (strcmp (ref, a1Seq)==0) {		//ref-alt on a2
			fprintf (des, "\t%s\t%s", ref, a2Seq);
			a1Gt = '0';
			a2Gt = '1';
		    }
		    if (strcmp (ref, a2Seq)==0) {		//ref-alt on a1
			fprintf (des, "\t%s\t%s", ref, a1Seq);
			a1Gt = '1';
			a2Gt = '0';
		    }
		    if (strcmp(ref, a1Seq)&&strcmp(ref, a2Seq)) {		//alt-alt 
			fprintf (des, "\t%s\t%s,%s", ref, a1Seq,a2Seq);
			a1Gt = '1';
			a2Gt = '2';
		    }
		  }//no no-calls
		 } 
		}
		} //not a hap call

		fprintf (des, "\t%d", QUAL);   
		fprintf (des, "\tPASS");   
		fprintf (des, "\tNS=1");   
		if (f_tRC) fprintf (des, ";DP=%d", tRC);   

		if (f_xRef) fprintf (des, ";DB");   

		if (f_tRC) strcpy (FORMAT, "GT:GQ:DP:HDP:HQ");
		else strcpy (FORMAT, "GT:GQ:HQ");

		if (strcmp (zg, "hap")==0) {					//haploid genome
		    if (f_tRC) snprintf (F1, 2000, "%c:%d:%d:%d,%d:%d,%d", a1Gt, QUAL, tRC, a1RC, 0, a1Sc, 0);
		    else snprintf (F1, 2000, "%c:%d:%d,%d", a1Gt, QUAL, a1Sc, 0);
		}
		else {								//diploid genome

	  	if (f_a1HL||f_a2HL) {                                            //phased
		    strcat (FORMAT, ":PS");
		    if (f_tRC) snprintf (F1, 2000, "%c|%c:%d:%d:%d,%d:%d,%d", a1Gt, a2Gt, QUAL, tRC, a1RC, a2RC, a1Sc, a2Sc);
		    else snprintf (F1, 2000, "%c|%c:%d:%d,%d", a1Gt, a2Gt, QUAL, a1Sc, a2Sc);
		    snprintf (Fa, 20000, ":%s,%s", a1HL, a2HL);
		    strcat (F1, Fa);
		}
		else {								//not phased
		    if (f_tRC) snprintf (F1, 2000, "%c/%c:%d:%d:%d,%d:%d,%d", a1Gt, a2Gt, QUAL, tRC, a1RC, a2RC, a1Sc, a2Sc);
		    else snprintf (F1, 2000, "%c/%c:%d:%d,%d", a1Gt, a2Gt, QUAL, a1Sc, a2Sc);
		}

		if (f_a1G||f_a2G) {
		    strcat (FORMAT, ":GENE");
		    kk = 0;
		    while (a1G[kk]!='\0') {
			if (a1G[kk]==':') a1G[kk]=',';
			kk ++;
		    }
		    kk = 0;
		    while (a2G[kk]!='\0') {
			if (a2G[kk]==':') a2G[kk]=',';
			kk ++;
		    }
		    snprintf (Fa, 20000, ":%s,%s", a1G, a2G);
		    strcat (F1, Fa);
		}
		}							//diploid regions

		if (f_mRNA) {
		    strcat (FORMAT, ":mRNA");
		    kk = 0;
		    while (mRNA[kk]!='\0') {
			if (mRNA[kk]==':') mRNA[kk]=',';
			kk ++;
		    }
		    snprintf (Fa, 20000, ":%s", mRNA);
		    strcat (F1, Fa);
		}

		if (f_rmsk) {
		    strcat (FORMAT, ":rmsk");
		    kk = 0;
		    while (rmsk[kk]!='\0') {
			if (rmsk[kk]==':') rmsk[kk]=',';
			kk ++;
		    }
		    snprintf (Fa, 20000, ":%s", rmsk);
		    strcat (F1, Fa);
		}

		if (f_segDup) {
		    strcat (FORMAT, ":segDup");
		    snprintf (Fa, 20000, ":%s", segDup);
		    strcat (F1, Fa);
		}
		if (f_rCov) {
		    strcat (FORMAT, ":rCov");
		    snprintf (Fa, 20000, ":%s", rCov);
		    strcat (F1, Fa);
		}
		if (f_cPd) {
		    strcat (FORMAT, ":cPd");
		    snprintf (Fa, 20000, ":%s", cPd);
		    strcat (F1, Fa);
		}

		fprintf (des, "\t%s\t%s\n", FORMAT, F1);

		pbg = bg;
	    }//ref != '=', variants presented
        }

        fclose (srV);
        fclose (srR);
        fclose (des);

return 0;
}


