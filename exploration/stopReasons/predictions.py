#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 07:44:19 2021

@author: olesyar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 12:15:38 2021

@author: olesyar
"""
# Libraries
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader, SequentialSampler
from transformers import BertModel
from transformers import BertTokenizer
from common_classes import BertClassifier
from common_classes import text_preprocessing
from common_classes import preprocessing_for_bert
from numpy import argmax
from common_classes import get_class
from common_classes import bert_predict
import csv
import torch.nn as nn
import logging
import datetime

import numpy as np
from pyspark.sql import SparkSession
import pyspark.sql.functions as F

CHEMBL_EVIDENCE_PATH = 'data/chembl-2021-08-23.json.gz'

model = torch.load('model/bert_trials.pth')
model.eval()

# Load the BERT tokenizer
tokenizer = BertTokenizer.from_pretrained('bert-base-uncased', do_lower_case=True)


logging.basicConfig(level=logging.ERROR)

# Columns of the studies file
names_studies = ['nct_id','nlm_download_date_description',
                 'study_first_submitted_date','results_first_submitted_date','disposition_first_submitted_date',
                 'last_update_submitted_date','study_first_submitted_qc_date','study_first_posted_date',
                 'study_first_posted_date_type','results_first_submitted_qc_date','results_first_posted_date',
                 'results_first_posted_date_type','disposition_first_submitted_qc_date',
                 'disposition_first_posted_date','disposition_first_posted_date_type',
                 'last_update_submitted_qc_date','last_update_posted_date','last_update_posted_date_type',
                 'start_month_year','start_date_type','start_date','verification_month_year',
                 'verification_date','completion_month_year','completion_date_type','completion_date',
                 'primary_completion_month_year','primary_completion_date_type','primary_completion_date',
                 'target_duration','study_type','acronym','baseline_population','brief_title','official_title',
                 'overall_status','last_known_status','phase','enrollment','enrollment_type','source',
                 'limitations_and_caveats','number_of_arms','number_of_groups','why_stopped','has_expanded_access',
                 'expanded_access_type_individual','expanded_access_type_intermediate',
                 'expanded_access_type_treatment','has_dmc','is_fda_regulated_drug','is_fda_regulated_device',
                 'is_unapproved_device','is_ppsd','is_us_export','biospec_retention','biospec_description',
                 'ipd_time_frame','ipd_access_criteria','ipd_url','plan_to_share_ipd','plan_to_share_ipd_description',
                 'created_at','updated_at']

def data_loader(study_stop_reason: np.ndarray):
    # Run `preprocessing_for_bert` on the column containing all the reasons to stop
    data_tokens, data_masks = preprocessing_for_bert(study_stop_reason)
    
    # Create the DataLoader for the data

    # Format tokens with attentions into Tensor to feed into the model
    data_input = TensorDataset(data_tokens, data_masks)
    data_sampler = SequentialSampler(data_input)
    data_dataloader = DataLoader(data_input, sampler=data_sampler, batch_size=32)
    print('The test set is ready')
    
    return data_dataloader

class Dataset(torch.utils.data.Dataset):

    def __init__(self, df):

        self.labels = [labels[label] for label in df['category']]
        self.texts = [tokenizer(text, 
                               padding='max_length', max_length = 512, truncation=True,
                                return_tensors="pt") for text in df['text']]

    def classes(self):
        return self.labels

    def __len__(self):
        return len(self.labels)

    def get_batch_labels(self, idx):
        # Fetch a batch of labels
        return np.array(self.labels[idx])

    def get_batch_texts(self, idx):
        # Fetch a batch of inputs
        return self.texts[idx]

    def __getitem__(self, idx):

        batch_texts = self.get_batch_texts(idx)
        batch_y = self.get_batch_labels(idx)

        return batch_texts, batch_y

# =============================================================================
# make predictions
# =============================================================================

def main():

    stopReasons = (
        spark.read.json(CHEMBL_EVIDENCE_PATH)

        # Extract a test set
        .sample(0.01)

        # Extract studies with their reasons to stop
        .withColumn('urls', F.explode('urls'))
        .filter(F.col('urls.niceName').contains('ClinicalTrials'))
        .withColumn('nct_id', F.element_at(F.split(F.col('urls.url'), '%22'), -2))
        .select('nct_id', 'studyStopReason')
        .filter(F.col('studyStopReason').isNotNull())
        .distinct()

        # Convert to Pandas DF
        #.toPandas()
    )
    logging.info('The test set is loaded.')

    # Preprocess study stop reasons
    study_stop_reason_array = (
        stopReasons.select('studyStopReason').toPandas()['studyStopReason'].values
    )

    probs = bert_predict(model, data_loader(study_stop_reason_array))
    csv_file1=open(f'output/stopped_predictions2-{datetime.today()}.tsv', "w")
    writer1 = csv.writer(csv_file1, delimiter='\t', lineterminator='\n')
    i=0
    not_stopped=stopReasons[stopReasons["why_stopped"].notnull()]    
    for ind,dat in not_stopped.iterrows():
        writer1.writerow([dat['why_stopped'].replace('\r~', ''),dat['phase'],dat['nct_id'],dat['start_date'], dat['overall_status'],dat['last_update_posted_date'],dat['completion_date'],get_class(argmax(probs[i]))])
        i=i+1
    csv_file2=open(f'output/notstopped_predictions2-{datetime.today()}.tsv', "w")
    writer2 = csv.writer(csv_file2, delimiter='\t')
    i=0
    stopped=stopReasons[stopReasons["why_stopped"].isnull()]
    for ind,dat in stopped.iterrows():
        writer2.writerow([dat['why_stopped'],dat['phase'],dat['nct_id'],dat['start_date'],dat['overall_status'],dat['last_update_posted_date'],dat['completion_date'],''])

if __name__ == '__main__':

    global spark
    spark = (
        SparkSession.builder
        .master('local[*]')
        .config("spark.driver.memory", "15g").appName('spark')
        .getOrCreate()
    )

