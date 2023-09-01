
import pandas as pd
import numpy as np
import torch
import torch.optim as optim
import pdb
from scipy.stats import beta

# def likelihood(X, W, Y, p, beta1, beta2, alpha):
#     logits = beta1 * X + beta2 * W
#     likelihoods = Y * torch.log(torch.sigmoid(logits)) + \
#                   (1 - Y) * torch.log(1 - torch.sigmoid(logits))
#     bernoulli_probs = torch.where(X == 1, p, 1 - p)  # p inference
#     
#     gaussian_probs = torch.exp(-0.5 * (W - (alpha * X)) ** 2) / \
#                      torch.sqrt(2 * torch.tensor(np.pi))
#     likelihoods += torch.log(gaussian_probs)
#     likelihoods += torch.log(bernoulli_probs)
#     
#     return torch.sum(likelihoods)

def likelihood(X1, X2, mu, a_alph):
    
    gaussian_probs = torch.exp(-0.5 * ((X1 - mu) / 40) ** 2) / \
                     torch.sqrt(40 * 2 * torch.tensor(np.pi))
    likelihoods = torch.log(gaussian_probs)
    # likelihoods += torch.log(bernoulli_probs)
    
    return torch.sum(likelihoods)

def optimize_parameters(dataset):
    X1 = torch.tensor(dataset['AverageMInFile'].values, dtype=torch.float32)
    X2 = torch.tensor(dataset['PercentTradesNeverDelq'].values, dtype=torch.float32)

    # initialize the parameters
    mu = torch.tensor(20, requires_grad=True, dtype=torch.float32)
    alpha0 = torch.tensor(0, requires_grad=True, dtype=torch.float32)
    # alpha1 = torch.tensor(0, requires_grad=True, dtype=torch.float32)
    # beta2 = torch.tensor(0, requires_grad=True, dtype=torch.float32)
    # alpha = torch.tensor(0, requires_grad=True, dtype=torch.float32)

    optimizer = optim.Adam([mu], lr=1)

    for _ in range(200):
        optimizer.zero_grad()
        loss = -likelihood(X1, mu) # neg-likelihood
        print('{:.3f}'.format(loss))
        print()
        loss.backward()
        optimizer.step()
        print('{:.3f}'.format(mu))
        print()

    return mu.item() # p.item(), alpha.item(), beta1.item(), beta2.item()

data = pd.read_csv('data/heloc_top6.csv', usecols=range(1, 7))

optimize_parameters(data)

