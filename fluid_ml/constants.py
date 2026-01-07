"""This module contains constants used throughout package."""
# Dictionary keys for the configuration files
PATHS = "paths"
TRAIN_FILES_PATH = "train_files_path"
SAVE_PATH = "save_path"
MODEL = "model"
TRAIN_RATIO = "train_ratio"
LEARNING_RATE = "learning_rate"
WEIGHT_DECAY = "weight_decay"
USE_LR_SCHEDULER = "use_lr_scheduler"
USE_EARLY_STOPPING = "use_early_stopping"
EARLY_STOPPING_PATIENCE = "early_stopping_patience"
CRITERION = "criterion"
BATCH_SIZE = "batch_size"
EPOCHS = "epochs"
NUM_HIDDEN_LAYERS = "num_hidden_layers"
KERNEL_SIZE = "kernel_size"
IN_CHANNELS = "in_channels"
OUT_CHANNELS = "out_channels"
HIDDEN_CHANNELS = "hidden_channels"
OUTPUT_ACTIVATION = "output_activation"
USE_BIAS = "use_bias"
PADDING_MODE = "padding_mode"
INPUTS = "inputs"
LABELS = "labels"
U = "u"  # noqa: WPS111
V = "v"  # noqa: WPS111
RANDOM_SPLIT = "random_split"
RANDOM_SEED = "random_seed"
MIN = "min"
MAX = "max"


# Paths
RESOURCES = "resources"
INPUTS_FILE_NAME = "inputs.pt"
BUILD = "build"
TRAIN = "train"
MODELS = "models"
DEFAULT_MODEL_NAME = "model.pt"
CONFIG_FILE_NAME = "config.json"
INIT_MODEL_PY = "mymodel.py"
MIN_MAX_YAML = "min_max.yaml"

# Plotting constants
FIG_HEIGHT = 5
FIG_WIDTH = 5
TITLE_FONT_SIZE = 20
TITLE_Y_POSITION = 0.95
IMG_SIZE = 21
N_CHANNELS = 2
COLORMAP = "coolwarm"

# Default training parameters
DEFAULT_TRAIN_RATIO = 0.8
DEFAULT_LR = 5e-5
DEFAULT_EARLY_STOPPING_PATIENCE = 20
DEFAULT_BATCH_SIZE = 32
DEFAULT_EPOCHS = 5000
DEFAULT_NUM_HIDDEN_LAYERS = 5
DEFAULT_KERNEL_SIZE = 7
DEFAULT_IN_CHANNELS = 1
DEFAULT_OUT_CHANNELS = 2
DEFAULT_HIDDEN_CHANNELS = 16
DEFAULT_OUTPUT_ACTIVATION = None
DEFAULT_USE_BIAS = True
DEFAULT_PADDING_MODE = "zeros"
DEFAULT_RANDOM_SEED = 0

# Normalization constants
EPSILON = 1e-9
