#!/usr/local/bin/python3.7 -u
import lib_process_data as lpd
import numpy
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--test')
    args = parser.parse_args()

    return


if __name__ == '__main__':
    lpd.print_script_call()
    main()
    lpd.save_script_call()
