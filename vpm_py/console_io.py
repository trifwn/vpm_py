import os

# Color codes
color_dict = {
    "red": "\033[91m",
    "green": "\033[92m",
    "blue": "\033[94m",
    "yellow": "\033[93m",
    "magenta": "\033[95m",
    "cyan": "\033[96m",
}
end_color = "\033[00m"

def print_blue(text, rank = 0):
    if rank == 0:
        print(f"\033[94m{text}\033[00m")

def print_green(text, rank = 0):
    if rank == 0:
        print(f"\033[92m{text}\033[00m")

def print_red(text, rank = 0):
    if rank == 0:
        print(f"\033[91m{text}\033[00m")

def print_IMPORTANT(
    text, 
    rank = 0,
    color_divider = "yellow",
    color_text = "red"
):
    # Get the color from the dictionary
    color_1 = color_dict[color_divider]
    color_2 = color_dict[color_text]
    if rank == 0:
        # try:
        #     _, columns = os.popen('stty size', 'r').read().split()
        # except:
        columns = 150
        print(f"{color_1}{'-'*int(columns)}{end_color}")
        # Print on the middle of the screen
        # If the line is multiline containing '\n', split it and print each line
        for line in text.split('\n'):
            print(f"{color_2}{line.center(int(columns))}{end_color}")
        print(f"{color_1}{'-'*int(columns)}{end_color}")

def print_divider(rank = 0):
    "Fill the screen with a divider"
    # Get the terminal size
    rows, columns = os.popen('stty size', 'r').read().split()
    if rank == 0:
        print(f"\033[93m{'-'*int(columns)}\033[00m")