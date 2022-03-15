from __future__ import print_function
from pyspark.sql import SparkSession
from pyspark import SparkContext, SparkConf
import numpy as np
import time
from operator import add
from huckel import huckel_2

#spark = SparkSession.builder.appName('Tugas Akhir').getOrCreate()
#sc = spark.sparkContext

conf = SparkConf().setAppName('Tugas Akhir').setMaster("spark://sparks.us-east1-c.c.sparkproject-319117.internal:7077")
# conf = SparkConf().setAppName('Tugas Akhir')
sc = SparkContext(conf=conf)
# spark = SparkSession.builder.appName('Tugas Akhir').getOrCreate()
# sc = spark.sparkContext

def main():
    t0 = time.perf_conter()
    x = np.arange(6,100,2)
    sc.parallelize(x).map(huckel_2).reduce(add)
    tn = time.perf_counter()
    print(f'Waktu komputasi = {tn-t0} s')
    sc.stop()

if __name__=='__main__':
    main()
