import os, subprocess

import cerberus
from cerberus.communicators import InetSocketCommunicator, UnixSocketCommunicator
from cerberus.solvers import Solver, CodeInput

# Logger options, maximum verbosity
cerberus.LOG.set_verbosity(3)

# Create the solver and initialize it. Decompose the case also if requested
parallel_run = True
decompose = True
n_proc = 8

path = os.environ["FOAMJOB_EXE_PATH"]
case_files = ["0", "constant", "system"]
subprocess.run(["./Allinit"])
if parallel_run:
    if decompose:
        subprocess.run(["foamDictionary", "-entry", "numberOfSubdomains", "-set", str(n_proc), "system/decomposeParDict"])
        subprocess.run(["decomposePar", "-force", "-allRegions"])
    proc_dirs = [f"processor{x}" for x in range(n_proc)]
    foam_input = CodeInput(case_files + proc_dirs, main_input_idx=None)
    foam = Solver("OF", path, ["-parallel", "krakenMultiRun"])
else:
    foam_input = CodeInput(case_files, main_input_idx=None)
    foam = Solver("OF", path, ["krakenMultiRun"])


foam.input = foam_input
foam.initialize(socket_comm_class=UnixSocketCommunicator)

# Output fields to communicate on each time step
output_fields = ["foam_of_density_fluid",
                 "foam_of_temperature_fluid",
                 "foam_of_density_solid",
                 "foam_of_temperature_solid"]

# Put 0.1 W power to each cell of the solid region
of_solid_power = foam.get_input_field("foam_if_power_solid")
of_solid_power.set_and_send(0.1)

# Time loop
nTimeSteps = 1000
dTime = 0.001
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

    # Get updated field values for the output fields
    for field_name in output_fields:
        tr = foam.get_output_field(field_name)
        tr.communicate()
