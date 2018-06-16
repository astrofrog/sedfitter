from __future__ import print_function, division

import os
import sys
import shutil

from astropy.table import Table

from ..six.moves import input

__all__ = ['create_dir', 'delete_dir', 'delete_file', 'read_table']


def create_dir(dir_name):
    delete_dir(dir_name)
    os.system("mkdir " + dir_name)


def delete_dir(dir_name):

    if os.path.exists(dir_name):
        reply = input("Delete directory " + dir_name + "? [y/[n]] ")
        if reply == 'y':
            shutil.rmtree(dir_name)
        else:
            print("Aborting...")
            sys.exit()
        print("")


def delete_file(file_name):

    if os.path.exists(file_name):
        reply = input("Delete file " + file_name + "? [y/[n]] ")
        if reply == 'y':
            os.remove(file_name)
        else:
            print("Aborting...")
            sys.exit()
        print("")


def read_table(file_name):
    try:
        return Table.read(file_name, format='fits', character_as_bytes=False)
    except TypeError:  # older versions of Astropy
        return Table.read(file_name, format='fits')
