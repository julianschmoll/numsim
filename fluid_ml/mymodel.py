# as in example jupyter notebook. this should probably be saved alongside the model
# maybe we never want this file here but only alongside the model to test
import torch

def init_my_model(in_channels:int=1,out_channels:int=2):
    model = torch.nn.Sequential()
    return model
