#!"i:\csc3002 - computer science project\csc3002\csc3002\environment\scripts\python.exe"
# EASY-INSTALL-ENTRY-SCRIPT: 'pathogenie','console_scripts','pathogenie-gui'
__requires__ = 'pathogenie'
import re
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(
        load_entry_point('pathogenie', 'console_scripts', 'pathogenie-gui')()
    )
