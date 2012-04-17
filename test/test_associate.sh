#!/bin/bash
if [[ $# == 0 ]]
then
  echo "./x.sh [table] [vcf_filename.gz] [project_filename.zip] [Rscript_filename.R] [random_num]"
  exit
fi

###
# source function
###
xpose () {
  linenum=`wc -l < $1`
  perl -lne '@x=split/\t/; foreach (0..$#x) {
  print "$_\t$.\t$x[$_]"}' $1 | \
  sort -k1n -k2n | \
  perl -ne 'chomp; @x=split/\t/,$_,3; print $x[2];
  print $.%'$linenum'?"\t":"\n"' 
}

unzip -q $3 -d tmp
cd tmp
logfile=`ls *.proj | sed 's/.proj/.log/g'`
###
# prepare input
###

# get genotype information from vtools tped
vtools export $1 --format tped --samples 1 --style numeric |\
  cut -f 1,4- > tped.geno.tmp
# add sample names && generate header line
sed "1i CHROM\tPOS\t`grep $2 $logfile | cut -f 2 | xargs echo`" tped.geno.tmp |\
  awk '{if(NR ==1) print "sample_name\t" $0; else print $1"_"$2"\t" $0}' |\
  sed 's/\s/\t/g' | cut -f1,4- > geno.tmp
# transpose to R input
xpose geno.tmp > tgeno.tmp

if [[ $5 == 0 ]]
then

	###
	# SNV analysis 
	###
	
	# By R
	join phenotype.tmp tgeno.tmp | sed 's/\s/\t/g' | Rscript $4 testSnv |\
	  # trim the output format 
	  sed 's/^X//g' | ./gw_round.py | grep -v -P "\tNA\t|pval" | sort > R.res
	
	# By vtools 
	vtools associate $1 BMI --covariate aff sex -m "LNBT --alternative 2" -j8 | \
	    awk '{if(NR ==1) print "ID\t" $0; else print $1"_"$2"\t" $0}' | sed 's/#chr/chr/g' | \
	    cut -f 1,5,6,7 | grep -v -P "\tNAN|pval" | sort > vtools.res

else

	###
	# association for groupby
	###
	
	# output allele information
	vtools output $1 chr pos ref alt --header |\
	# add random group
	  ./fakecols.py $5 > vtools.group.tmp 
	  # add new group in vtools
	vtools update $1 --from_file vtools.group.tmp --format cache/vtools_randcol.fmt --var_info grpby
	# vtools association by group
	vtools associate $1 BMI --covariate aff sex -m "LNBT --alternative 2" -g grpby -j8 |\
	    cut -f 1,3,4,5 | grep -v -P "\tNA|pval" | sort -g > vtools.res
        # permutation based approach for empirical p-values
        vtools associate $1 BMI --covariate aff sex -m "LNBT --alternative 2 -p 5000 --permute_by X --adaptive 0.00001" -g grpby -j8 |\
            grep -v -P "\tNA|pval" | cut -f 1,3 | sort -g > vtools2.res	
        # vtools association by variable threshold method
        vtools associate $1 BMI --covariate aff sex -m "VariableThresholdsQt --alternative 2 -p 5000 --permute_by X --adaptive 0.00001" -g grpby -j8 |\
            grep -v -P "\tNA|pval" | cut -f 1,3 | awk '{if ($2>=0) $2=$2; else $2=0-$2; print $0}' | sed 's/\s/\t/g' | sort -g > vtools3.res
        #compare t values
        cut -f1,2 vtools.res > vtools1.res
        diff vtools1.res vtools2.res
        join vtools3.res vtools2.res | sed 's/\s/\t/g' | awk '{if ($2>$3) exit; else if($2=$3) exit; else if($2<$3) print $1}'
	###
	# group analysis by R
	###
	
	awk '{if(NR ==1) print "sample_name\t" $0; else print $1"_"$2"\t" $0}' vtools.group.tmp | cut -f 1,6 > group.tmp
	join group.tmp geno.tmp | sed 's/\s/\t/g' | Rscript $4 scoreRegion |\
        awk '{if(NR ==1) print "sample_name\t" $0; else print "" $0}' > R.group.tmp
	join phenotype.tmp R.group.tmp | sed 's/\s/\t/g' | Rscript $4 testRegion | sed 's/^X//g' |\
        ./gw_round.py | grep -v -P "\tNA|pval" | sort -g > R.res
fi

#compare R result to vtools result
diff R.res vtools.res

###
# clear up history
###
cd ..
rm -r tmp
