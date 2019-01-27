y1file=$1
x1file=$2
outdir=$3

/net/gs/vol1/home/alexhu/proj/siren/bin/C++/ratest_batch -p $x1file -q $y1file -a both -o $outdir -i 500 -s 0 -u -z 0.01 # The -u is to make it unsquared!! take heed!

