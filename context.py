import os
import getpass

if getpass.getuser() == "kadu":
    home_dir = "/home/kadu/Dropbox/splus-halpha"
elif getpass.getuser() == "55119":
    home_dir = "/c/Users/55119/dropbox/splus-halpha"

data_dir = os.path.join(home_dir, "data")
if not os.path.exists(data_dir):
    os.mkdir(data_dir)

bands = ['U', 'F378', 'F395', 'F410', 'F430', 'G', 'F515', 'R', 'F660', 'I',
         'F861', 'Z']
<<<<<<< HEAD
=======

tables_dir = os.path.join(home_dir, "tables")
>>>>>>> e0e612189c11c14b4a0399fad07aedaff0e53a0e
