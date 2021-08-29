## ReVal-SNV

setwd('codes')

if (!dir.exists('../results')) dir.create('../results')

#### train models
source('get_VM_est.R')

#### make calls
source('get_VM_calls.R')
