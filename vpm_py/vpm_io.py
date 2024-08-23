import os

def print_blue(text, rank = 0):
    if rank == 0:
        print(f"\033[94m{text}\033[00m")

def print_green(text, rank = 0):
    if rank == 0:
        print(f"\033[92m{text}\033[00m")

def print_red(text, rank = 0):
    if rank == 0:
        print(f"\033[91m{text}\033[00m")

def print_IMPORTANT(text, rank = 0):
    if rank == 0:
        print(f"\033[93m{'-'*100}\033[00m")
        print(f"\033[91m{text}\033[00m")
        print(f"\033[93m{'-'*100}\033[00m")
