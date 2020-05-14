# -*- coding: UTF-8 -*-

"""
base class and functions used in training, prediction, plotting
"""
import numpy as np
import pandas as pd
from PIL import Image
from pathlib import Path
import torch
import torch.nn as nn
from torch.utils.data import Dataset
from torchvision import transforms, models

def split_df_to2(df, n):
    '''
    args:
        df: pandas dataframe
        n (int): the number of rows for part one
    return:
        df1: dataframe of part1
        df2: dataframe of part2
        
    '''
    df1 = df.sample(n)
    df2 = df.loc[df.index[~df.index.isin(df1.index)], :]
    return df1, df2

def split_val_train(df, n_val, n_fold):
    '''
    args:
        df: the dataframe for training and validataion for one category
        n_val (int): number of validation data in each fold
    yield fold, index of training data, index of validation data
    '''
    all_idx = df.index
    if n_val*n_fold > all_idx.shape[0]:
        print('not enough data for %s fold cross-validation!'%n_fold)
        return
    val_idx_st, val_idx_ed = 0, n_val
    for fold in range(1, n_fold+1):
        val_idx = all_idx[val_idx_st: val_idx_ed]  
        train_idx = all_idx[~all_idx.isin(val_idx)]
        val_idx_st += n_val
        val_idx_ed += n_val
        yield fold, train_idx, val_idx

class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""
    def __init__(self, mn_prefix, patience=20, verbose=True, delta=0):
        """
        Args:
            mn_prefix (str): the prefix of the saved model name.
            patience (int): How long to wait after last time validation loss improved.
                            Default: 20
            verbose (bool): If True, prints a message for each validation loss improvement. 
                            Default: True
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            Default: 0
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta
        self.mn_prefix = mn_prefix

    def __call__(self, val_loss, model):
        score = -val_loss

        if self.best_score is None: # the first epoch
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        '''Saves model when validation loss decrease.'''
        if self.verbose:
            print(f'Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
        torch.save(model.state_dict(), '%s.pt'%self.mn_prefix)
        self.val_loss_min = val_loss

class LeafcountingDataset(Dataset):
    """leaf counting dataset."""

    def __init__(self, csv_file, root_dir, transform=None):
        """
        Args:
            csv_file (string): Path to the comma separated csv file without header. The 1st column is image file name and the 2nd column is the annotation/label. 
            root_dir (string): Directory with all the images.
        """
        self.csv = pd.read_csv(csv_file, header=None)
        self.root_dir = Path(root_dir)
        self.transform = transform

    def __len__(self):
        return len(self.csv)

    def __getitem__(self, idx):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        img_name = self.csv.iloc[idx, 0]
        image = Image.open(self.root_dir/img_name)
        if len(image.getbands()) == 4:
            image = image.convert('RGB')
        label = self.csv.iloc[idx, 1].astype('float32').reshape(-1,)

        if self.transform:
            image = self.transform(image)

        return image, label, img_name

image_transforms = {
    # Train uses data augmentation
    'train':
    transforms.Compose([
        transforms.RandomRotation(degrees=15),
        transforms.ColorJitter(),
        transforms.RandomHorizontalFlip(),
        transforms.Resize(size=(224, 224)),
        transforms.ToTensor(),
        transforms.Normalize([0.485, 0.456, 0.406],
                             [0.229, 0.224, 0.225])  # Imagenet standards
    ]),
    # Validation does not use augmentation
    'valid':
    transforms.Compose([
        transforms.Resize(size=(224, 224)),
        transforms.ToTensor(),
        transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])
    ]),
}

def train_model_regression(model, dataloaders, criterion, optimizer, model_name_prefix, patience=10, num_epochs=10, is_inception=False):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device == 'cuda:0':
        model.cuda()
    valid_loss_history = []
    train_loss_history = []
    
    early_stopping = EarlyStopping(model_name_prefix, patience=patience, verbose=True)
    for epoch in range(num_epochs):
        print('Epoch {}/{}'.format(epoch, num_epochs - 1))
        print('-' * 10)

        # Each epoch has a training and validation phase
        for phase in ['train', 'valid']:
            if phase == 'train':
                model.train()  # Set model to training mode
            else:
                model.eval()   # Set model to evaluate mode

            running_loss = 0.0
            for inputs, labels, img_fn in dataloaders[phase]:
                inputs = inputs.to(device)
                labels = labels.to(device)
                
                # zero the parameter gradients
                optimizer.zero_grad()

                # forward
                # track history if only in train
                with torch.set_grad_enabled(phase == 'train'):
                    if is_inception and phase == 'train':
                        # From https://discuss.pytorch.org/t/how-to-optimize-inception-model-with-auxiliary-classifiers/7958
                        outputs, aux_outputs = model(inputs)
                        loss1 = criterion(outputs, labels)
                        loss2 = criterion(aux_outputs, labels)
                        loss = loss1 + 0.4*loss2
                    else:
                        outputs = model(inputs)
                        loss = criterion(outputs, labels)
                    # backward + optimize only if in training phase
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()
                # statistics
                running_loss += loss.item() * inputs.size(0)
            epoch_loss = running_loss / len(dataloaders[phase].dataset)
            print('{} Loss: {:.4f}'.format(phase, epoch_loss))
            if phase == 'train':
                train_loss_history.append(epoch_loss)
            if phase == 'valid':
                valid_loss_history.append(epoch_loss)
                early_stopping(epoch_loss, model)
        if early_stopping.early_stop:
            print('Early stopping')
            break
    print('Best val loss: {:4f}'.format(early_stopping.val_loss_min))

    # load best model weights
    model.load_state_dict(torch.load('%s.pt'%model_name_prefix))
    return model, train_loss_history, valid_loss_history

def set_parameter_requires_grad(model, feature_extracting):
    if feature_extracting:
        for param in model.parameters():
            param.requires_grad = False

def initialize_model(model_name='vgg16', num_classes=1, feature_extract=True, use_pretrained=True):
    # Initialize these variables which will be set in this if statement. Each of these
    #   variables is model specific.
    model_ft = None
    input_size = 0

    if model_name == "googlenet":
        """ googlenet
        """
        model_ft = models.googlenet(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = 224

    elif model_name == "resnet152":
        """ Resnet152
        """
        model_ft = models.resnet152(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = 224

    elif model_name == "vgg16":
        """ VGG16
        """
        model_ft = models.vgg16(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs,num_classes)
        input_size = 224

    elif model_name == "resnet18":
        """ Resnet18
        """
        model_ft = models.resnet18(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs, num_classes)
        input_size = 224

    elif model_name == "alexnet":
        """ Alexnet
        """
        model_ft = models.alexnet(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs,num_classes)
        input_size = 224

    elif model_name == "vgg11_bn":
        """ VGG11_bn
        """
        model_ft = models.vgg11_bn(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier[6].in_features
        model_ft.classifier[6] = nn.Linear(num_ftrs,num_classes)
        input_size = 224

    elif model_name == "squeezenet":
        """ Squeezenet
        """
        model_ft = models.squeezenet1_0(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        model_ft.classifier[1] = nn.Conv2d(512, num_classes, kernel_size=(1,1), stride=(1,1))
        model_ft.num_classes = num_classes
        input_size = 224

    elif model_name == "densenet":
        """ Densenet
        """
        model_ft = models.densenet121(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        num_ftrs = model_ft.classifier.in_features
        model_ft.classifier = nn.Linear(num_ftrs, num_classes)
        input_size = 224

    elif model_name == "inception":
        """ Inception v3
        Be careful, expects (299,299) sized images and has auxiliary output
        """
        model_ft = models.inception_v3(pretrained=use_pretrained)
        set_parameter_requires_grad(model_ft, feature_extract)
        # Handle the auxilary net
        num_ftrs = model_ft.AuxLogits.fc.in_features
        model_ft.AuxLogits.fc = nn.Linear(num_ftrs, num_classes)
        # Handle the primary net
        num_ftrs = model_ft.fc.in_features
        model_ft.fc = nn.Linear(num_ftrs,num_classes)
        input_size = 299

    else:
        print("Invalid model name, exiting...")
        exit()

    return model_ft, input_size