__author__ = 'vasilev_is'
import pydoc
import os

files = [f for f in os.listdir('.') if os.path.isfile(f)]
files = [f.replace('.py','') for f in files if 'Fianora' in f and 'Estimator_Decorator' not in f]

pydoc.writedoc('Fianora_Derivative')

# for f in files:
#     #pydoc.help(f)
#
#     pydoc.writedoc(f)

import shutil

#s = [f for f in os.listdir('.') if os.path.isfile(f) and 'html' in f]

[print (f) for f in files]

# for f in files:
#     shutil.move (f,  'pydocs/')
#
#
