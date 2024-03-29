import os
import getpass
from sys import argv
import matplotlib.pyplot as plt


home_dir= []
if getpass.getuser() == "kadu":
    home_dir = "/home/kadu/Dropbox/splus-halpha"
elif getpass.getuser() == "amori":
    home_dir = "c:/Users/amori/Dropbox/splus-halpha (1)"

data_dir = os.path.join(home_dir, "data")
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
    print(data_dir)

bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
         'F861', 'Z']

wave_eff = {"F378": 3773.4, "F395": 3940.8, "F410": 4895.4,
            "F430": 4292.5, "F515": 5133.5, "F660": 6614.0, "F861": 8608.1,
            "G": 4647.8, "I": 7683.8, "R": 6266.6, "U": 3536.5,
            "Z": 8679.5}

tables_dir = os.path.join(home_dir, "tables")

# Matplotlib settings
plt.style.context("seaborn-paper")
plt.rcParams["text.usetex"] = False
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['font.serif'] = 'Computer Modern'
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True

SMALL_SIZE = 7
MEDIUM_SIZE = 8
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# set tick width
width = 0.5
majsize = 4
minsize = 2
plt.rcParams['xtick.major.size'] = majsize
plt.rcParams['xtick.major.width'] = width
plt.rcParams['xtick.minor.size'] = minsize
plt.rcParams['xtick.minor.width'] = width
plt.rcParams['ytick.major.size'] = majsize
plt.rcParams['ytick.major.width'] = width
plt.rcParams['ytick.minor.size'] = minsize
plt.rcParams['ytick.minor.width'] = width
plt.rcParams['axes.linewidth'] = width