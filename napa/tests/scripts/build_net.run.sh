napa=$1
config=$2

mkdir -p log
logfile = log/build.$config.oe

echo "Logging in:" 
echo $logfile

python $napa/run_napa.py \
    -r build \
    -c $config >& $logfile &
