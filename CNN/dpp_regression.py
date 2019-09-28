"""
train neural network to detect whether plant flowers or not
"""
import warnings
warnings.filterwarnings('ignore',category=FutureWarning)
from glob import glob
import numpy as np
import pickle
import deepplantphenomics as dpp
from pathlib import Path
import os
import sys

def train(train_dir, label_fn, model_dir, epoch, lr):
    """
    train_dir: the directory where your training images located
    label_fn: the file name of labels under train_dir. Just specify the file name don't inclde the path. 
    model_dir: the name of you model. Model results will be save to this dir
    epoch: specify the epoch. Based on dpp document suggest 100 for plant stress and 500 for counting.
    lr: specify learnning rate. 0.0001 used in dpp leaf counting example
    """
    model_dir_path = Path(model_dir)
    if not model_dir_path.exists():
        model_dir_path.mkdir()
    tensorboard_dir_path = model_dir_path/'tensorboard'
    img_dir = Path(train_dir)

    model = dpp.RegressionModel(debug=True, save_checkpoints=True, report_rate=150, tensorboard_dir=str(tensorboard_dir_path), save_dir=str(model_dir_path))
    model.set_batch_size(45)
    model.set_number_of_threads(10)
    model.set_image_dimensions(418, 283, 3)
    model.set_resize_images(True)

    model.set_num_regression_outputs(1)
    model.set_test_split(0.0)
    model.set_validation_split(0.0)
    model.set_learning_rate(float(lr))
    model.set_weight_initializer('xavier')
    model.set_maximum_training_epochs(int(epoch))

    # Augmentation options
    model.set_augmentation_brightness_and_contrast(True)
    model.set_augmentation_flip_horizontal(True)
    model.set_augmentation_flip_vertical(True)
    #model.set_augmentation_crop(True)
    
    # Load labels and images
    model.load_multiple_labels_from_csv(img_dir/label_fn, id_column=0)
    model.load_images_with_ids_from_directory(img_dir)

    # Define a model architecture
    model.add_input_layer()
    model.add_convolutional_layer(filter_dimension=[5, 5, 3, 32], stride_length=1, activation_function='tanh')
    model.add_pooling_layer(kernel_size=3, stride_length=2)
    model.add_convolutional_layer(filter_dimension=[5, 5, 32, 64], stride_length=1, activation_function='tanh')
    model.add_pooling_layer(kernel_size=3, stride_length=2)
    model.add_convolutional_layer(filter_dimension=[3, 3, 64, 64], stride_length=1, activation_function='tanh')
    model.add_pooling_layer(kernel_size=3, stride_length=2)
    model.add_convolutional_layer(filter_dimension=[3, 3, 64, 64], stride_length=1, activation_function='tanh')
    model.add_pooling_layer(kernel_size=3, stride_length=2)
    model.add_output_layer()
    # Begin training the model
    model.begin_training()

if len(sys.argv)==6:
    train(*sys.argv[1:])
else:
    print('train_dir', 'label_fn', 'model_dir', 'epoch', 'lr')
