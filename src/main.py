# This si the main script

import argparse
import NFS
import cProfile

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', required=True, help='number to factor')
    
    args = parser.parse_args()
    
    n = args.n

    NFS.NFS(n)