#!"i:\csc3002 - computer science project\csc3002\csc3002\environment\scripts\python.exe"
# EASY-INSTALL-ENTRY-SCRIPT: 'ete3==3.1.2','console_scripts','ete3'
__requires__ = 'ete3==3.1.2'
import re
import sys
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(
        load_entry_point('ete3==3.1.2', 'console_scripts', 'ete3')()
    )
