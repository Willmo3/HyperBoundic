#!/usr/bin/env python3
import os
import subprocess
import time
from sys import stdout

# check if the script is run from the correct directory
if not os.path.isfile('run_sanity_tests.py'):
    print("Error: This script must be run from the directory containing 'run_sanity_tests.py'.")
    exit(1)

# check if the executable is placed in out
if not os.path.isfile('out/PDEapprox'):
    print("Error: The executable 'out/PDEapprox' was not found. Please build the project first.")
    exit(1)

# generate files
subprocess.run(['./out/PDEapprox', '-w'], check=True)

# Now, run sanity tests for all permutations of files
domains = ['real', 'interval', 'affine', 'mixed']
fluxes = ['cubic', 'burgers', 'lwr', 'buckley_leverett']
solvers = ['lax_friedrichs', 'leapfrog']

for domain in domains:
    for flux in fluxes:
        conds = f'simulations/{flux}_{domain}_conds.json'
        for solver in solvers:
            print(f'Running test: Domain={domain}, Flux={flux}, Solver={solver}')
            print()
            stdout.flush() # ensure output is printed in order if redirected to file
            config = f'simulations/{domain}_{flux}_{solver}_config.json'

            cmd = f'./out/PDEapprox -c {config} -s {conds}'

            time_before = time.perf_counter()
            subprocess.run(['./out/PDEapprox', '-c', config, '-s', conds], check=True)
            time_after = time.perf_counter()
            print()
            print(f'Test completed in {time_after - time_before:.4f} seconds.\n')
            stdout.flush() # ensure output is printed in order if redirected to file