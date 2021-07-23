from huckel import huckel_2
import numpy as np
from numpy import linalg as la
import time
from scipy import linalg as la
import threading

def main():
    start = time.time()
    X = np.arange(6,101,2)
    for x in X:
        t1 = threading.Thread(target=huckel_2, args=(x,))
        t1.start()
        t1.join()
    stop = time.time()
    print('Waktu komputasi = {} s'.format(stop-start))

if __name__ == "__main__":
    main()