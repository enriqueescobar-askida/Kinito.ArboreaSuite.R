
# seriously, the input file are in en_CA.ISO8859 and they might contains ponctuation and accents

export LANG=en_CA.ISO8859

mkdir -p qc
mkdir -p normalize
mkdir -p results
mkdir -p control_analysis

export DISPLAY=:0.0

R --vanilla < script.R > script.R.STDOUT 2> script.R.STDERR

