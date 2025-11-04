#!/usr/bin/env python3
import os
import subprocess
from sys import stdout

# check if the script is run from the correct directory
if not os.path.isfile('run_sanity_tests.py'):
    print("Error: This script must be run from the directory containing 'run_sanity_tests.py'.")
    exit(1)

# check if the executable is placed in out
if not os.path.isfile('out/PDEapprox'):
    print("Error: The executable 'out/executable' was not found. Please build the project first.")
    exit(1)

# generate files
subprocess.run(['./out/PDEapprox', '-w'], check=True)

# Now, run sanity tests for all permutations of files
domains = ['real', 'interval', 'affine', 'mixed']
fluxes = ['cubic', 'burgers', 'lwr', 'buckley_leverett']
solvers = ['lax_friedrichs', 'leapfrog']

for domain in domains:
    conds = f'simulations/{domain}_conds.json'
    for flux in fluxes:
        for solver in solvers:
            print(f'Running test: Domain={domain}, Flux={flux}, Solver={solver}')
            stdout.flush() # ensure output is printed in order if redirected to file
            config = f'simulations/{domain}_{flux}_{solver}_config.json'

            cmd = f'./out/PDEapprox -c {config} -s {conds}'
            subprocess.run(['./out/PDEapprox', '-c', config, '-s', conds], check=True)
            print()
            stdout.flush() # ensure output is printed in order if redirected to file