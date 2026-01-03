"""This module contains constants used throughout package."""
# Dictionary keys for the configuration files
PATHS_KEY = "paths"
TRAIN_FILES_PATH_KEY = "train_files_path"
BASE_SAVE_PATH_KEY = "save_path"
MODEL_SAVE_PATH_KEY = "model"
TRAIN_RATIO_KEY = "train_ratio"
LEARNING_RATE_KEY = "lr"
WEIGHT_DECAY_KEY = "weight_decay"
USE_LR_SCHEDULER_KEY = "use_lr_scheduler"
USE_EARLY_STOPPING_KEY = "use_early_stopping"
EARLY_STOPPING_PATIENCE_KEY = "early_stopping_patience"
CRITERION_KEY = "criterion"
BATCH_SIZE_KEY = "batch_size"
EPOCHS_KEY = "epochs"
NUM_HIDDEN_LAYERS_KEY = "num_hidden_layers"
KERNEL_SIZE_KEY = "kernel_size"
IN_CHANNELS_KEY = "in_channels"
OUT_CHANNELS_KEY = "out_channels"
HIDDEN_CHANNELS_KEY = "hidden_channels"
OUTPUT_ACTIVATION_KEY = "output_activation"
USE_BIAS_KEY = "use_bias"
PADDING_MODE_KEY = "padding_mode"
INPUTS_KEY = "inputs"
LABELS_KEY = "labels"
U_CHANNEL_KEY = "u"
V_CHANNEL_KEY = "v"
RANDOM_SPLIT_KEY = "random_split"
MIN = "min"
MAX = "max"


# Paths
RESOURCE_FOLDER_NAME = "resources"
INPUTS_FILE_NAME = "inputs.pt"
BUILD_PATH = "build"
TRAIN_FOLDER_NAME = "train"
MODELS_SAVE_PATH = "models"
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

# Normalization constants
EPSILON = 1e-9
