# !/bin/bash
if [[ $# == 0 ]]
then
  echo "./x.sh [table] [vcf_filename.gz] [project_filename.zip] [random_num]"
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

xpose_std () {
  data=`< /dev/stdin`
  linenum=`echo $data | wc -l`
  echo $data | \
  perl -lne '@x=split/\t/; foreach (0..$#x) {
  print "$_\t$.\t$x[$_]"}' | \
  sort -k1n -k2n | \
  perl -ne 'chomp; @x=split/\t/,$_,3; print $x[2];
  print $.%'$linenum'?"\t":"\n"' 
}

zcatvcf () {
  zcat $1 | grep -v '^##' | grep -v '0|2'| grep -v '2|0' | cut -f 1,2,10- | sed 's/:[0-9]*:[A-Z]*//g' | sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/0\/1/1/g' | sed 's/1|1/2/g' | awk '{if(NR ==1) print "ID\t" $0; else print $1"_"$2"\t" $0}' | sed 's/#//g' | sed 's/\s/\t/g'
}

unzip -q $3 -d tmp
cd tmp

###
# prepare input
###

# get genotype information from vtools tped
vtools export $1 --format tped --samples "sample_name is not null" --style numeric |\
  cut -f 1,4- > tped.geno.tmp
# add sample names && generate header line
sed "1i CHROM\tPOS\t`grep $2 unitest.log | cut -f 2 | xpose_std`" tped.geno.tmp |\
  awk '{if(NR ==1) print "sample_name\t" $0; else print $1"_"$2"\t" $0}' |\
  sed 's/\s/\t/g' | cut -f1,4- > geno.tmp
# transpose to R input
xpose geno.tmp > tgeno.tmp

if [[ $4 == 0 ]]
then

	###
	# SNV analysis 
	###
	
	# By R
	join phenotype.tmp tgeno.tmp | sed 's/\s/\t/g' | Rscript test_associate.R testSnv |\
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
	  ./fakecols.py $4 > vtools.group.tmp 
	  # add new group in vtools
	vtools update $1 --from_file vtools.group.tmp --format cache/vtools_randcol.fmt --var_info grpby
	# vtools association by group
	vtools associate $1 BMI --covariate aff sex -m "LNBT --alternative 2" -g grpby -j8 |\
	    cut -f 1,3,4,5 | grep -v -P "\tNA|pval" | sort > vtools.res
	
	###
	# group analysis by R
	###
	
	awk '{if(NR ==1) print "sample_name\t" $0; else print $1"_"$2"\t" $0}' vtools.group.tmp | cut -f 1,6 > group.tmp
	join group.tmp geno.tmp | sed 's/\s/\t/g' | Rscript test_associate.R scoreRegion |\
        awk '{if(NR ==1) print "sample_name\t" $0; else print "" $0}' > R.group.tmp
	join phenotype.tmp R.group.tmp | sed 's/\s/\t/g' | Rscript test_associate.R testRegion | sed 's/^X//g' |\
        ./gw_round.py | grep -v -P "\tNA|pval" | sort > R.res
fi

#compare R result to vtools result
diff R.res vtools.res

###
# clear up history
###
cd ..
rm -r tmp
