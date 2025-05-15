#!/bin/bash --norc
#
function die {
    echo "$*" 1>&2
    exit 1
}

type perl > /dev/null || \
    die "$0: *** FATAL, this script requires program: perl"
pversion=$(perl -e 'print "$]\n"')
echo "Perl: using version $pversion"

# check installed Perl modules
perl_modules="Getopt::Long File::Basename MIME::Base64 \
              XML::Parser XML::LibXML XML::LibXSLT"
for pmod in $perl_modules
do
    pversion=$(perl -e "use $pmod; print \$$pmod::VERSION")
    if [ ${#pversion} -eq 0 ]; then
	die "$0: *** FATAL, perl module $pmod not installed"
    else
	echo "$pmod: using version $pversion"
    fi
done
echo -e "Good, all required modules are available on this system\n"

# Check validity of installation
[ ! -x ./hpf_eef_transform.pl ] && chmod +x ./hpf_eef_transform.pl
[ ! -d Parser_out_Ref ] && \
    die "$0: *** FATAL, directory ./Parser_out_Ref is missing"
[ ! -d Parser_out ] && \
    die "$0: *** FATAL, directory ./Parser_out is missing"
rm -f Parser_out/*

# perform transformation of Level 1b & 2 products
L1B_PROD_IDS="EGG_NOM_1b SST_NOM_1b STR_VC2_1b STR_VC3_1b"
L2_PROD_IDS="EGG_NOM_2_ EGG_TRF_2_ EGM_GOC_2_ SST_AUX_2_ SST_PSO_2_"
for PROD_ID in $L1B_PROD_IDS $L2_PROD_IDS 
do
    files=$(ls EEF/GO_CONS_${PROD_ID}_*)
    echo "Processing productID: $PROD_ID this may take a while..."
    if [ -x /usr/bin/time ]; then
	export TIME="Elapsed real time %E (%P), Memory used (max) %M KB\n"
	/usr/bin/time ./hpf_eef_transform.pl --output Parser_out $files
    else
	export TIMEFORMAT=$'Elapsed real time %E (%P%%)\n'
	time ./hpf_eef_transform.pl --output Parser_out $files
    fi
done

# check transformation of Level 1b & 2 products
declare -i nr_failed=0
files_ref=$(ls Parser_out_Ref/*)
for fl_ref in $files_ref
do 
    fl_out=Parser_out/$(basename $fl_ref)
    if [ ! -r $fl_out ]; then 
	echo "$0: *** FATAL, missing product: $(basename $fl_ref)"
	let nr_failed+=1 
    elif ! cmp $fl_ref $fl_out; then 
	echo "$0: *** FATAL, corrupted product: $(basename $fl_ref)" 
	let nr_failed+=1
    fi
done
[ $nr_failed -eq 0 ] && \
    echo -e "Good, all products have been sucessfully generated\n"

exit 0
