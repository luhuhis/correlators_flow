#!/usr/bin/env python3
import lib_process_data as lpd
import numpy
import argparse



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test')
    args = parser.parse_args()
    return args


def main():

    args = parse_args()

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
