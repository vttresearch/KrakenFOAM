import os, subprocess

import cerberus
from cerberus.solvers import Solver, CodeInput

# Logger options, maximum verbosity
cerberus.LOG.set_verbosity(3)

# Create the solver and initialize it. Decompose the case also if requested
parallel_run = True
decompose = True
n_proc = 4

path = os.environ["FOAMJOB_EXE_PATH"]
case_files = ["0", "constant", "system", "all_cell_centres"]
if parallel_run:
    if decompose:
        subprocess.run(["foamDictionary", "-entry", "numberOfSubdomains", "-set", str(n_proc), "system/decomposeParDict"])
        subprocess.run(["decomposePar", "-force"])
    proc_dirs = [f"processor{x}" for x in range(n_proc)]
    foam_input = CodeInput(case_files + proc_dirs, main_input_idx=None)
    foam = Solver("OF", path, ["-parallel", "krakenRun"])
else:
    foam_input = CodeInput(case_files, main_input_idx=None)
    foam = Solver("OF", path, ["krakenRun"])

foam.input = foam_input
foam.initialize()

# Get fields
of_power = foam.get_input_field("foam_if_power_core")
of_density = foam.get_output_field("foam_of_density_core")
of_temperature = foam.get_output_field("foam_of_temperature_core")

# Set power of 1e5 W in each cell
of_power.set_and_send(1.0e+05)

# Time loop
nTimeSteps = 3000
dTime = 0.005
startTime = 0.0

for tidx in range(nTimeSteps):

    t1 = (tidx+1)*dTime + startTime
    t0 = (tidx)*dTime + startTime

    # Set the time step
    if tidx == 0:
        foam.set_current_time_interval(t0, t1)
    else:
        foam.advance_to_time_interval(t0, t1)

    t0, t1 = foam.return_current_time_interval();
    print(f"Time interval = {t0}, {t1}")

    # Solve the time step
    foam.solve()

# Get temperature field from OpenFOAM and print the values to terminal
of_temperature.communicate()

print(f"Temperature field:\n")
print(of_temperature.value_vec)
