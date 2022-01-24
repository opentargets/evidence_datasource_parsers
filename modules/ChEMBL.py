#!/usr/bin/env python
"""This module adds the category of why a clinical trial has stopped early to the ChEMBL evidence."""

import argparse
import logging
import sys
from typing import Optional

from pyspark import SparkFiles
from pyspark.conf import SparkConf
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame
import pyspark.sql.functions as F

from common.evidence import detect_spark_memory_limit, write_evidence_strings

