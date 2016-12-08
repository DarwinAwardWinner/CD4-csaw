#!/usr/bin/env python

import dis
import os.path
import sys

from tempfile import TemporaryFile
from subprocess import check_output, Popen, PIPE
from importlib import import_module
python_files = check_output(''' find . -maxdepth 2 -name '*.py' -not -name list-python-modules.py ''', shell=True)\
               .decode(sys.getdefaultencoding()).strip().split('\n')
snakefiles = check_output(''' find . -maxdepth 2 -name '*Snakefile' ''', shell=True)\
               .decode(sys.getdefaultencoding()).strip().split('\n')

with TemporaryFile(mode='w+t') as tempf:
    for pf in python_files:
        tempf.write(open(pf, mode='rt').read())
        tempf.write('\n')
    for sf in snakefiles:
        p = Popen(['snakemake', '--print-compilation', '--snakefile', sf], stdout=PIPE)
        tempf.write(p.stdout.read().decode(sys.getdefaultencoding()))
        tempf.write('\n')
    tempf.flush()
    tempf.seek(0)
    fulltext = tempf.read()
instructions = filter(lambda x: x.opname == 'IMPORT_NAME', dis.get_instructions(fulltext))
imports = { x.argval for x in instructions }
# print('\n'.join(sorted(imports)))

import_files = {}
for libname in imports:
    lib = import_module(libname)
    try:
        import_files[libname] = lib.__file__
    except AttributeError:
        pass
print('\n'.join(''.join([k, ": ", repr(v)]) for k,v in import_files.items()))
