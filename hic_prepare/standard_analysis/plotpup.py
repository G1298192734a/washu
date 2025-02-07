#!/public/frasergen/PUB/software/python/Python-3.9.12/bin/python3
# -*- coding: utf-8 -*-
import re
import sys
from coolpuppy.plotpuppy_CLI import main
if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
