#!/bin/bash
set -e

GEN_CONFIG=$(readlink -f "./cfg/ml/train.txt")
TRAIN_CONFIG=$(readlink -f "./cfg/ml/config.json")

usage() {
    echo "Usage: $0 [-g generation_config] [-t training_config]"
    echo ""
    echo "Options:"
    echo "  -g    Path to the data generation config (Default: ./cfg/ml/train.txt)"
    echo "  -t    Path to the ML training config (Default: ./cfg/ml/config.json)"
    echo "  -h    Display this help message"
    exit 1
}

while getopts "g:t:h" opt; do
  case $opt in
    g) GEN_CONFIG=$(readlink -f "$OPTARG") ;;
    t) TRAIN_CONFIG=$(readlink -f "$OPTARG") ;;
    h) usage ;;
    \?) usage ;;
  esac
done

echo "----------------------------------------"
echo "SIMULATION CONFIG: $GEN_CONFIG"
echo "ML CONFIG:         $TRAIN_CONFIG"
echo "----------------------------------------"

echo "----------------------------------------"
echo "Building and Running Simulations..."
echo "----------------------------------------"

rm -rf build
mkdir -p build && cd build
cmake ..
make install
mpirun numsim_parallel "$GEN_CONFIG"
cd ..

echo "----------------------------------------"
echo "Training ML Model..."
echo "----------------------------------------"

python ./fluid_ml/main.py --config "$TRAIN_CONFIG"
