from huckel import huckel_2
import time
import multiprocessing as mp
import numpy as np

def main():
    start = time.time()
    X = np.arange(6,101,2)
    with mp.Pool(5) as p:
        res = p.map(huckel_2,X)
        print(res[0][3])
    stop = time.time()
    print('Waktu komputasi = {} s'.format(stop-start))

if __name__ == "__main__":
    main()
