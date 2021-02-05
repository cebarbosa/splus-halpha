import os
import getpass

if getpass.getuser() == "kadu":
    home_dir = "/home/kadu/Dropbox/splus-halpha"
elif getpass.getuser() == "55119":
    home_dir = "C:/Users/55119/Dropbox/splus-halpha (1)"

data_dir = os.path.join(home_dir, "data")
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
         'F861', 'Z']

tables_dir = os.path.join(home_dir, "tables")